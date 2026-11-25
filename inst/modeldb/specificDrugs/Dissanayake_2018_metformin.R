Dissanayake_2018_metformin <- function() {
  description <- "One-compartment oral population PK model for immediate-release metformin in metformin-naive adults with type 2 diabetes and chronic kidney disease stage 4, parameterised by the elimination rate constant, apparent clearance and absorption rate constant, with correlated additive-normal between-subject variability on CL/F and ka (nlme fit; no covariates) (Dissanayake 2018)"
  reference <- "Dissanayake AM, Wheldon MC, Hood CJ. Pharmacokinetics of metformin in patients with chronic kidney disease stage 4 and metformin-naive type 2 diabetes. Pharmacol Res Perspect. 2018;6(5):e00424. doi:10.1002/prp2.424"
  vignette <- "Dissanayake_2018_metformin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The random effects are additive-normal on the natural parameter scale: the
  # paper fitted Cl/F and ka 'as patient-specific random intercepts' in R nlme
  # and reports their SDs in L/h and 1/h (Table 2). They therefore act on the
  # linear typical values cl and ka, not on the log-scale intercepts lcl_int / lka_int.
  paper_specific_etas <- c("etacl", "etaka")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "metformin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "metformin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "30-75 years (eligibility)",
    age_median = "66.0 years (IQR 6.54)",
    weight_range = "77.2-149.8 kg (Table A1)",
    sex_female_pct = 22.2,
    disease_state = "Metformin-naive type 2 diabetes (diagnosed >= 2 years, HbA1c 6-11%) with chronic kidney disease stage 4 (eGFR 15-29 mL/min/1.73 m^2); not on renal replacement therapy",
    dose_range = "Immediate-release metformin hydrochloride 250, 500 or 1000 mg orally once daily (6 patients per dose cohort), taken fasting before breakfast",
    regions = "New Zealand",
    renal_function = "eGFR median 21.0 mL/min/1.73 m^2 (IQR 8.0); Cockcroft-Gault creatinine clearance on adjusted body weight 18-49 mL/min, median 30.2 mL/min (IQR 8.6) (Results; Table A1)",
    bmi = "Median 38.0 kg/m^2 (IQR 9.87); all patients >= 23 kg/m^2",
    sampling = "Single-dose profile at 0, 2, 4, 6, 8 and 24 h after the first dose (102 concentrations, 5-6 per patient), plus fasting troughs 24 h post-dose at later visits up to Day 32. The PK model was fitted to the single-dose profiles only. Serum metformin by HPLC, LOQ about 20 ug/L.",
    notes = "Phase I open-label consecutive-cohort dose-escalation study (Dissanayake 2018, Methods 'Subjects and protocol'). 4 of 18 patients female (2 in the 250 mg group, 2 in the 1000 mg group). Baseline demographics summarised in the Results text and Table A1."
  )

  ini({
    # Equation (1): c(t) = D*K*ka / ((Cl/F)*(ka - K)) * (exp(-K*t) - exp(-ka*t)),
    # estimated parameters K, ka and Cl/F; V/F = (Cl/F)/K is derived.
    lka_int <- log(0.65); label("Absorption rate constant ka (1/h)") # Table 2, ka fixed effect = 0.65 1/h (95% CI 0.45-0.86)
    lcl_int <- log(29.6); label("Apparent oral clearance CL/F (L/h)") # Table 2, Cl/F fixed effect = 29.6 L/h (95% CI 23.8-35.4)
    lkel <- log(0.119); label("Elimination rate constant K (1/h)") # Table 2, K fixed effect = 0.119 1/h (95% CI 0.11-0.13)

    # Additive-normal random effects on the natural scale (Table 2 'Random effects'):
    #   sd(Cl/F) = 11.56 L/h -> variance 133.6336
    #   sd(ka)   = 0.395 1/h -> variance 0.156025
    #   cor      = 0.785     -> covariance 0.785 * 11.56 * 0.395 = 3.584467
    etacl + etaka ~ c(133.6336, 3.584467, 0.156025) # Table 2, sd(Cl/F) = 11.56 L/h, sd(ka) = 0.395 1/h, cor(Cl/F,ka) = 0.785

    addSd <- 0.095; label("Additive residual error (mg/L)") # Table 2, 'Resid. SE' = 0.095 mg/L (95% CI 0.081-0.113)
  })
  model({
    ka <- exp(lka_int) + etaka
    cl <- exp(lcl_int) + etacl
    kel <- exp(lkel)
    # K carries no random effect, so V/F follows CL/F: vc = cl / kel
    # (typical value 29.6 / 0.119 = 248.7 L; Table 2 prints V/F = 249.7 L).
    vc <- cl / kel

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, vc in L: central / vc is mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
