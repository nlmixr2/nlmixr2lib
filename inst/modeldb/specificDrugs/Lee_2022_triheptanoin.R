Lee_2022_triheptanoin <- function() {
  description <- paste(
    "One-compartment population PK model for heptanoate, the main circulating",
    "metabolite of oral triheptanoin, in healthy adults and in pediatric and",
    "adult patients with long-chain fatty acid oxidation disorders (LC-FAOD)",
    "(Lee 2022). Each triheptanoin dose is split between two parallel",
    "first-order absorption compartments: depot (typical 56.2 percent, ka1",
    "0.425 1/h) and depot2 (43.8 percent, ka2 0.507 1/h, lag 3.74 h).",
    "Apparent clearance and volume scale allometrically with time-varying body",
    "weight (reference 58 kg); clearance is 19 percent lower in LC-FAOD",
    "patients than in healthy subjects and falls with time since the first",
    "dose by up to 45 percent (T50 86.7 h). Doses are entered as umol",
    "triheptanoin; each mole releases three moles of heptanoate."
  )
  reference <- paste(
    "Lee SK, Gosselin NH, Jomphe C, McKeever K, Putnam W. Population",
    "Pharmacokinetics of Heptanoate in Healthy Subjects and Patients With",
    "Long-Chain Fatty Acid Oxidation Disorders Treated With Triheptanoin.",
    "Clin Pharmacol Drug Dev. 2022;11(11):1264-1272. doi:10.1002/cpdd.1145.",
    "Final estimates from Table 3; bootstrap results and simulated exposures",
    "from Supplementary Tables S3 and S4.",
    sep = " "
  )
  vignette <- "Lee_2022_triheptanoin"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # An oral triheptanoin dose is given as TWO dose records of the full dose
  # (umol triheptanoin), one to depot and one to depot2; f(depot) and
  # f(depot2) then route F1 and 1 - F1 of it (Figure 1). Doses are in umol of
  # triheptanoin (MW 428.61 g/mol, Methods); the factor 3 in f() converts to
  # umol of heptanoate, as Table 3 footnote a applies to CL/F and V/F.
  dosing <- c("depot", "depot2")

  compartmentData <- list(
    depot = list(analyte = "heptanoate", units = "umol", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "heptanoate", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "heptanoate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying actual body weight at dosing (Table 3 footnote c). Allometric scaling of CL/F and V/F with estimated exponents, reference 58 kg (median WT of the PK population).",
      source_name = "WT"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with LC-FAOD)",
      notes = "Table 3 reports separate typical CL/F0 values for LC-FAOD subjects (471 L/h) and healthy subjects (584 L/h); encoded as a multiplicative ratio 584/471 applied when DIS_HEALTHY = 1. Retained in the final model for clinical relevance although not statistically significant (Results, Model Development).",
      source_name = "disease population (healthy vs LC-FAOD)"
    )
  )

  covariatesDataExcluded <- list(
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units = "mL/min/1.73 m^2 (age < 12 y, Schwartz) or mL/min (age >= 12 y, Cockcroft-Gault)",
      type = "continuous",
      notes = "Trend with the CL/F random effect (Figure S1); tested stepwise but not retained (Results, Model Development)."
    ),
    TBILI_BASE = list(
      description = "Baseline total bilirubin",
      units = "mg/dL (as reported; the register's canonical unit is umol/L)",
      type = "continuous",
      notes = "Trend with the CL/F random effect (Figure S1); tested stepwise but not retained (Results, Model Development)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 43,
    n_studies = 3,
    age_range = "0.87-62.1 years (LC-FAOD pediatric 0.87-17.9 y; LC-FAOD adult 19.2-62.1 y)",
    age_mean = "healthy 39.1 y; LC-FAOD pediatric 7.41 y; LC-FAOD adult 33.6 y",
    weight_range = "8.23-121 kg",
    weight_median = "58 kg (overall PK population)",
    sex_female_pct = 46.5,
    race_ethnicity = c(White = 81.4, Black = 9.3, Asian = 2.3, NativeHawaiianPacificIslander = 2.3, Other = 4.7),
    disease_state = "13 healthy adults and 30 patients with long-chain fatty acid oxidation disorders (23 pediatric, 7 adult)",
    dose_range = "Triheptanoin oil 1.25 or 1.5 g/kg/day orally in healthy subjects; 25-35 percent of daily caloric intake in at least four daily doses with food or by gastrostomy tube in LC-FAOD",
    regions = "not reported",
    n_observations = 562,
    notes = paste(
      "Table 1 (continuous) and Table 2 (categorical) baseline demographics;",
      "Supplementary Table S1 lists the three studies (phase 1 healthy",
      "crossover, phase 2 LC-FAOD study, long-term extension PK sub-study)",
      "and Supplementary Table S2 the 562 included heptanoate concentrations",
      "(LLOQ 1.0 umol/L). Two subjects contributed to two studies under",
      "separate IDs, so the 43 analysed IDs correspond to 41 people."
    )
  )

  ini({
    # Absorption (Table 3). F1 is held on the logit scale (Table 3 footnote b:
    # F1 = exp(tvF1 + nF) / (1 + exp(tvF1 + nF))). It is encoded here as the
    # fraction routed to the delayed compartment, frac = 1 - F1, following the
    # registered dual-absorption idiom; logit(1 - F1) = -logit(F1), and the
    # symmetric eta keeps the same variance, so the two forms are identical.
    lka <- log(0.425); label("First-order absorption rate constant from depot, Ka1 (1/h)") # Table 3: Ka1 = 0.425 1/h
    lka2 <- log(0.507); label("First-order absorption rate constant from depot2, Ka2 (1/h)") # Table 3: Ka2 = 0.507 1/h
    ltlag2 <- log(3.74); label("Absorption lag time of depot2, Lag2 (h)") # Table 3: Lag2 = 3.74 h
    logitfrac <- log((1 - 0.562) / 0.562); label("Logit of the fraction of the dose absorbed from depot2, 1 - F1 (unitless)") # Table 3: F1 = 0.562

    # Disposition (Table 3; reference weight 58 kg). Footnote c:
    # CL/F = CL/F0 * (WT/58)^Expo_CL * (1 - MAX * TIME / (T50 + TIME)).
    lcl <- log(471); label("Apparent clearance CL/F0 for a 58 kg LC-FAOD patient at time 0 (L/h)") # Table 3: CL/F0 = 471 L/h (LC-FAOD subjects)
    e_healthy_cl <- 584 / 471; label("Ratio of healthy-subject to LC-FAOD CL/F0, applied as ratio^DIS_HEALTHY (unitless)") # Table 3: CL/F0 = 584 L/h (healthy subjects) vs 471 L/h (LC-FAOD)
    e_wt_cl <- 1.07; label("Allometric exponent of body weight on CL/F, Expo_CL (unitless)") # Table 3: Expo_CL = 1.07
    lcl_time_max <- log(0.450); label("Maximum fractional decrease of CL/F with time since first dose, MAX (unitless)") # Table 3: MAX = 0.450
    lcl_t50 <- log(86.7); label("Time to half of the maximum CL/F decrease, T50 (h)") # Table 3: T50 = 86.7 (h; footnote c 'TIME is time in hours')
    lvc <- log(5.56); label("Apparent volume of distribution V/F0 for a 58 kg subject (L)") # Table 3: V/F0 = 5.56 L
    e_wt_vc <- 1.13; label("Allometric exponent of body weight on V/F, Expo_V (unitless)") # Table 3: Expo_V = 1.13

    # Between-subject variability. Table 3 footnote b: the F1 entry is the
    # variance on the logit scale; the others are %BSV = sqrt(exp(omega^2) - 1).
    etalogitfrac ~ 1.09 # Table 3: BSV F1 = 1.09 (logit-scale variance; shrinkage 46.7%)
    etalka ~ log(1 + 1.34^2) # Table 3: BSV Ka1 = 134% (shrinkage 45.0%)
    etalcl ~ log(1 + 0.601^2) # Table 3: BSV CL/F = 60.1% (shrinkage 26.7%)

    # Residual error: log-additive (Methods), reported as the SD on the log scale.
    expSd <- 0.764; label("Log-additive residual error SD (log scale)") # Table 3: Log-additive error = 0.764
  })

  model({
    # Individual absorption parameters
    ka <- exp(lka + etalka)
    ka2 <- exp(lka2)
    tlag2 <- exp(ltlag2)
    frac <- expit(logitfrac + etalogitfrac)

    # Disposition (Table 3 footnotes c and d). `t` is time since the first
    # triheptanoin dose in hours.
    cl_time_max <- exp(lcl_time_max)
    cl_t50 <- exp(lcl_t50)
    cl <- exp(lcl + etalcl) * e_healthy_cl^DIS_HEALTHY * (WT / 58)^e_wt_cl *
      (1 - cl_time_max * t / (cl_t50 + t))
    vc <- exp(lvc) * (WT / 58)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central

    # One mole of triheptanoin yields three moles of heptanoate (Table 3
    # footnote a); the dose split F1 / (1 - F1) is Figure 1.
    f(depot) <- 3 * (1 - frac)
    f(depot2) <- 3 * frac
    alag(depot2) <- tlag2

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
