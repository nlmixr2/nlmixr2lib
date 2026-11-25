Aruldhas_2019_rifampicin <- function() {
  description <- paste(
    "One-compartment population PK model for oral rifampicin in Indian",
    "children aged 2-16 years with pulmonary or lymph-node tuberculosis",
    "(Aruldhas 2019), sampled at steady state (fully auto-induced) in the",
    "intensive phase of daily or thrice-weekly therapy. Absorption through",
    "nine fixed transit compartments with the absorption rate set equal to",
    "the transit rate (ka = ktr = 10 / MTT). Fixed allometric scaling of",
    "CL/F (exponent 0.75) and V/F (exponent 1) on body weight centred at",
    "19.4 kg. Relative bioavailability fixed to 1 with between-subject",
    "variability; BSV also on V/F and MTT, none on CL/F; additive error on",
    "log concentrations. Companion isoniazid model:",
    "modellib('Aruldhas_2019_isoniazid').",
    sep = " "
  )
  reference <- paste(
    "Aruldhas BW, Hoglund RM, Ranjalkar J, Tarning J, Mathew SK,",
    "Verghese VP, Bose A, Mathew BS. Optimization of dosing regimens of",
    "isoniazid and rifampicin in children with tuberculosis in India.",
    "Br J Clin Pharmacol. 2019;85(3):644-654. doi:10.1111/bcp.13846.",
    sep = " "
  )
  vignette <- "Aruldhas_2019_tuberculosis"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit9 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric function on CL/F (exponent 0.75) and V/F",
        "(exponent 1) (Methods, 'Population pharmacokinetic modelling';",
        "Results: estimating the exponent did not improve the fit). Centred",
        "at 19.4 kg, the 'typical child' of the Table 4 footnote. Table 2",
        "reports median (IQR) weight 19.5 (13.7-33.7) kg in the 41 enrolled",
        "children.",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  # Screened in the stepwise covariate search and not retained ('No other
  # covariates could be retained in the stepwise covariate approach').
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods); not retained."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Evaluated as an allometric size descriptor alternative to body",
        "weight (Methods); body weight was retained.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 39L,
    n_studies = 1L,
    n_observations = 284L,
    age_range = "2-16 years; median 7 (IQR 3.5-13) in the 41 enrolled children (Table 2)",
    weight_range = "median 19.5 kg (IQR 13.7-33.7) in the 41 enrolled children (Table 2)",
    sex_female_pct = 29.3,
    race_ethnicity = "Indian (single centre, Vellore, Tamil Nadu)",
    disease_state = paste(
      "Newly diagnosed pulmonary (n = 36) or lymph-node (n = 5)",
      "tuberculosis, haemoglobin >= 10 g/dL; relapse, multidrug-resistant",
      "tuberculosis and HIV co-infection excluded. Sampled during the last",
      "2 weeks of the 8-week intensive phase (steady state, after about",
      "60 days of rifampicin, so auto-induction is complete), after an",
      "overnight fast.",
      sep = " "
    ),
    dose_range = paste(
      "Rifampicin with isoniazid, pyrazinamide and ethambutol by the",
      "Revised National Tuberculosis Control Programme weight bands",
      "(Table 1, Regimen 2: 75-300 mg per dose), given daily (WHO,",
      "n = 14) or thrice weekly (RNTCP directly observed therapy,",
      "n = 27).",
      sep = " "
    ),
    regions = "India (Christian Medical College, Vellore, Tamil Nadu)",
    notes = paste(
      "41 children enrolled; 2 were excluded from the rifampicin analysis",
      "for an implausible ~2 h absorption delay, leaving 284",
      "concentrations from 39 children (16% below the 0.04 mg/L LLOQ,",
      "imputed as LLOQ/2). Samples pre-dose and at 0.5, 1, 1.5, 2, 2.5, 4",
      "and 6 h. Auto-induction could not be estimated from steady-state",
      "data, so the model describes the induced state only.",
      sep = " "
    )
  )

  ini({
    # Table 4 population estimates are 'for a typical child weighing 19.4 kg
    # with tuberculosis' (Table 4 footnote).

    # Absorption: nine fixed transit compartments with ka = ktr (Results:
    # 'kA and kTR were set to be equal, since estimating them separately did
    # not improve the model fit'). The mean transit time spans the 9 + 1 = 10
    # first-order transfers from depot to central, so ktr = 10 / MTT. The
    # vignette confirms this reproduces the published median Tmax of 1.59 h.
    lmtt <- log(0.932) ; label("Mean absorption transit time MTT (h)") # Table 4 row 'MTT (h)' = 0.932 (RSE 9.92%)

    lcl <- log(8.11) ; label("Apparent clearance CL/F (L/h, 19.4-kg child)") # Table 4 row 'CL/F (l h-1)' = 8.11 (RSE 10.9%)
    lvc <- log(44.7) ; label("Apparent volume of distribution V/F (L, 19.4-kg child)") # Table 4 row 'Vc/F (l)' = 44.7 (RSE 14.6%)

    # 'Relative bioavailability (F) was fixed to unity for the population but
    # with an estimated interindividual variability' (Methods).
    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability F (unitless)") # Table 4 row 'F (%)' = 100 fix

    # 'Effect of both body weight and fat-free mass were evaluated with
    # allometric scaling on all clearances and volumes parameters with an
    # exponent of 0.75 and 1, respectively' (Methods); estimating the
    # exponent did not improve the fit (Results).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL/F (unitless)") # Methods, 'Population pharmacokinetic modelling': exponent 0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of body weight on V/F (unitless)") # Methods, 'Population pharmacokinetic modelling': exponent 1

    # Between-subject variability. Table 4 reports '%CV for BSV'; converted
    # with omega^2 = log(1 + CV^2), the exact log-normal relation. No BSV was
    # estimated on CL/F (Table 4 BSV column is blank for the CL row).
    etalvc ~ 0.1625 # Table 4 row 'Vc/F' BSV = 42% CV (RSE 12.6%); shrinkage 3.40%
    etalmtt ~ 0.2410 # Table 4 row 'MTT (h)' BSV = 52.2% CV (RSE 12.7%); shrinkage 5.06%
    etalfdepot ~ 0.3843 # Table 4 row 'F (%)' BSV = 68.0% CV (RSE 16.6%); shrinkage 10.71%

    # Residual error: 'additive error on the individually predicted
    # logarithmic concentrations' (Methods), i.e. lnorm() here. The Table 4
    # RUV of 0.271 is read as the VARIANCE sigma^2: its bootstrap 95% CI
    # (0.191-0.379) is twice as wide, relative to the estimate, as its
    # printed RSE of 8.96% implies -- the signature of a variance whose RSE
    # is reported on the SD scale. The BSV rows, which are on the CV scale,
    # show no such factor. See the vignette.
    expSd <- sqrt(0.271) ; label("Additive residual error on log concentration (SD)") # Table 4 row 'RUV' = 0.271 (RSE 8.96%), read as a variance
  })
  model({
    # Allometric size scaling centred at the 19.4-kg typical child.
    cl <- exp(lcl) * (WT / 19.4)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 19.4)^e_wt_vc
    mtt <- exp(lmtt + etalmtt)

    # Nine transit compartments; ka = ktr, so MTT covers 10 transfers.
    ktr <- 10 / mtt
    ka <- ktr
    kel <- cl / vc

    d / dt(depot) <- -ktr * depot
    d / dt(transit1) <- ktr * depot - ktr * transit1
    d / dt(transit2) <- ktr * transit1 - ktr * transit2
    d / dt(transit3) <- ktr * transit2 - ktr * transit3
    d / dt(transit4) <- ktr * transit3 - ktr * transit4
    d / dt(transit5) <- ktr * transit4 - ktr * transit5
    d / dt(transit6) <- ktr * transit5 - ktr * transit6
    d / dt(transit7) <- ktr * transit6 - ktr * transit7
    d / dt(transit8) <- ktr * transit7 - ktr * transit8
    d / dt(transit9) <- ktr * transit8 - ktr * transit9
    d / dt(central) <- ka * transit9 - kel * central

    f(depot) <- exp(lfdepot + etalfdepot)

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
