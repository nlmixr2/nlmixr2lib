Dissanayake_2018_metformin_crcl <- function() {
  description <- "One-compartment oral population PK model for immediate-release metformin in metformin-naive adults with type 2 diabetes and chronic kidney disease stage 4, with apparent clearance proportional to baseline Cockcroft-Gault creatinine clearance (CL/F = CrCl x a CL/F-to-CrCl ratio) and correlated additive-normal between-subject variability on that ratio and on ka (nlme fit; the paper's dose-recommendation model) (Dissanayake 2018)"
  reference <- "Dissanayake AM, Wheldon MC, Hood CJ. Pharmacokinetics of metformin in patients with chronic kidney disease stage 4 and metformin-naive type 2 diabetes. Pharmacol Res Perspect. 2018;6(5):e00424. doi:10.1002/prp2.424"
  vignette <- "Dissanayake_2018_metformin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The random effects are additive-normal on the natural parameter scale: the
  # paper fitted (Cl/F)' and ka as patient-specific random effects in R nlme
  # and reports their SDs on the parameter scale (Table A3). etae_crcl_cl
  # therefore acts on the linear ratio e_crcl_cl and etaka on the linear ka.
  paper_specific_etas <- c("etae_crcl_cl", "etaka")

  covariateData <- list(
    CRCL = list(
      description = "Baseline creatinine clearance by the modified Cockcroft-Gault equation, RAW and NOT BSA-normalised. The weight in the Cockcroft-Gault numerator is chosen by BMI band: total body weight if BMI < 18.5 kg/m^2, ideal body weight (Devine: 50 kg for men or 45.5 kg for women + 2.3 kg per inch of height over 60 inches) if BMI 18.5-22.9 kg/m^2, and adjusted body weight IBW + 0.4 x (TBW - IBW) if BMI >= 23 kg/m^2. Every study patient had BMI >= 23 kg/m^2, so the adjusted-body-weight column of Table A1 is the value used.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (baseline). The paper fitted the model with CrCl in L/h; model() converts the mL/min column with CRCL * 60 / 1000. Enters as a through-origin multiplier on CL/F with no centring: CL/F = CrCl (L/h) x 17.4. Cohort median 30.2 mL/min (IQR 8.6), range 17.96-48.55 mL/min (Results; Table A1). The authors judged the model's steady-state predictions unreliable below about 30-32 mL/min and gave dose recommendations only for 30-50 mL/min (Table A4).",
      source_name = "CrCl"
    )
  )

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
    notes = "Phase I open-label consecutive-cohort dose-escalation study (Dissanayake 2018, Methods 'Subjects and protocol'). 4 of 18 patients female. AIC/BIC preferred the companion model without CrCl (Dissanayake_2018_metformin; AIC -57.0 vs -48.9), but the authors used this model for the CrCl-specific dose recommendations of Table A4."
  )

  ini({
    # Equation (2): c(t) = D*K*ka / (CrCl*(Cl/F)'*(ka - K)) * (exp(-K*t) - exp(-ka*t)),
    # with CrCl in L/h and (Cl/F)' = (Cl/F) / CrCl.
    lka_int <- log(0.564); label("Absorption rate constant ka (1/h)") # Table A3, ka fixed effect = 0.564 1/h (95% CI 0.401-0.729)
    lkel <- log(0.128); label("Elimination rate constant K (1/h)") # Table A3, K fixed effect = 0.128 1/h (95% CI 0.113-0.143)
    e_crcl_cl <- 17.4; label("Ratio of CL/F to creatinine clearance, both in L/h (unitless)") # Table A3, Cl/F' fixed effect = 17.4 (95% CI 13.5-21.3); Results 'Pharmacokinetics and kidney function'

    # Additive-normal random effects on the natural scale (Table A3 'Random effects'):
    #   sd(Cl/F') = 7.81  -> variance 60.9961
    #   sd(ka)    = 0.31  -> variance 0.0961
    #   cor       = 0.79  -> covariance 0.79 * 7.81 * 0.31 = 1.912669
    etae_crcl_cl + etaka ~ c(60.9961, 1.912669, 0.0961) # Table A3, sd(Cl/F') = 7.81, sd(ka) = 0.31 1/h, cor(Cl/F',ka) = 0.79

    # Table A3 prints 'Resid. SE' = 0.99 mg/L with 95% CI (0.08, 0.12); the
    # point estimate lies outside its own interval and is ten-fold the
    # companion model's 0.095 mg/L, so the decimal point is misplaced and the
    # value is 0.099 mg/L.
    addSd <- 0.099; label("Additive residual error (mg/L)") # Table A3, 'Resid. SE' printed 0.99 with 95% CI 0.08-0.12; read as 0.099
  })
  model({
    crcl_lh <- CRCL * 60 / 1000

    ka <- exp(lka_int) + etaka
    cl <- crcl_lh * (e_crcl_cl + etae_crcl_cl)
    kel <- exp(lkel)
    # K carries no random effect, so V/F follows CL/F: vc = cl / kel.
    vc <- cl / kel

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, vc in L: central / vc is mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
