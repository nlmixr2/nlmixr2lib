Tanigawara_2012_garenoxacin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "garenoxacin in Japanese healthy volunteers and patients with",
    "community-acquired pneumonia or secondary infection of chronic",
    "respiratory disease (Tanigawara 2012 final model, phase I-III data).",
    "Apparent clearance is linear-additive in Cockcroft-Gault creatinine",
    "clearance and body weight; apparent volume is a power function of body",
    "weight plus an additive shift for men."
  )
  reference <- paste(
    "Tanigawara Y, Nozawa K, Tsuda H. Optimal dose finding of garenoxacin",
    "based on population pharmacokinetics/pharmacodynamics and Monte Carlo",
    "simulation. Eur J Clin Pharmacol. 2012;68(1):39-53.",
    "doi:10.1007/s00228-011-1095-3"
  )
  vignette <- "Tanigawara_2012_garenoxacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "garenoxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "garenoxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation on total body weight, raw mL/min (NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tanigawara 2012 Table 2 footnote b: male CLCR = (140 - age) x WT /",
        "(72 x SCr in mg/dL); female value times 0.85. Linear-additive effect",
        "on CL/F centred at 85 mL/min, the arithmetic mean of the pooled",
        "phase I-III analysis set (Table 2: 102.4 mL/min in the 104 phase I/II",
        "subjects and 71.1 mL/min in the 133 phase III patients). Observed",
        "range 19.3-199.2 mL/min."
      ),
      source_name = "CLCR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear-additive effect on CL/F centred at 55 kg and power effect on",
        "V/F normalised to 55 kg, the arithmetic mean of the pooled analysis",
        "set (Table 2: 58.5 kg phase I/II, 52.7 kg phase III). Observed range",
        "33-120 kg."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0) carries the additive V/F shift",
      notes = paste(
        "Tanigawara 2012 codes GEN = 1 for male and 0 for female (Table 3",
        "footnote a) and adds theta6 x GEN to V/F. Stored here as the",
        "canonical SEXF with GEN = 1 - SEXF, so the coefficient keeps its",
        "published sign: men have a V/F 6.98 L larger than women of the same",
        "weight."
      ),
      source_name = "GEN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 237L,
    n_studies = 5L,
    n_observations = 1742L,
    age_range = "19-94 years",
    weight_range = "33-120 kg",
    sex_female_pct = 30.8,
    race_ethnicity = c(Japanese = 100),
    disease_state = paste(
      "49 healthy volunteers (three phase I studies), 55 patients with",
      "community-acquired pneumonia (CAP) or secondary infection of chronic",
      "respiratory disease (SICRD) in phase II, and 133 patients with mild or",
      "moderate SICRD in phase III."
    ),
    dose_range = paste(
      "Phase I 100-600 mg single dose, 200 or 400 mg once daily for 7-14 days,",
      "and a 400 mg fed/fasted study; phase II 200 or 400 mg once daily for",
      "7-14 days; phase III 400 mg once daily for 10 days. All oral."
    ),
    regions = "Japan",
    renal_function = "Cockcroft-Gault CLCR 19.3-199.2 mL/min (phase I/II mean 102.4, phase III mean 71.1)",
    notes = paste(
      "Tanigawara 2012 Tables 1-2. Men 164 of 237 (77/104 in phase I/II and",
      "87/133 in phase III). Dense sampling in phase I, sparse (1-4 samples",
      "per subject) in phases II and III. Final model estimated with NONMEM",
      "VI FOCE with interaction."
    )
  )

  ini({
    # Tanigawara 2012 Table 4, 'End of phase III stage' column (final model,
    # NONMEM VI FOCE-I). CL/F = theta1 + theta2 x (CLCR - 85) + theta3 x
    # (WT - 55) in mL/min; V/F = theta4 x (WT/55)^theta5 + theta6 x GEN in L.
    # The clearance terms are kept in the published mL/min units and converted
    # to L/h in model().
    lka <- log(2.14); label("First-order absorption rate constant (1/h)") # Table 4 final model: ka = 2.14 +/- 0.223 1/h
    lcl <- log(64.3); label("Apparent clearance CL/F at CLCR 85 mL/min and WT 55 kg (mL/min)") # Table 4 final model: theta1 = 64.3 +/- 1.25 mL/min
    lvc <- log(53.9); label("Apparent volume V/F for a 55 kg woman (L)") # Table 4 final model: theta4 = 53.9 +/- 1.89 L

    e_crcl_cl <- 0.213; label("Additive change in CL/F per mL/min of (CRCL - 85) (mL/min per mL/min)") # Table 4 final model: theta2 = 0.213 +/- 0.0309
    e_wt_cl <- 0.276; label("Additive change in CL/F per kg of (WT - 55) (mL/min per kg)") # Table 4 final model: theta3 = 0.276 +/- 0.113
    e_wt_vc <- 0.739; label("Power exponent of (WT / 55) on V/F (unitless)") # Table 4 final model: theta5 = 0.739 +/- 0.104
    e_sex_vc <- 6.98; label("Additive V/F shift for men, applied as (1 - SEXF) (L)") # Table 4 final model: theta6 = 6.98 +/- 1.88 L, GEN = 1 male

    # Table 4 final model omega^2 on the log scale (P_j = P x exp(eta_j), Eq. 2):
    # CL/F 0.0403 (CV 20.1%), V/F 0.0249 (CV 15.8%), covariance 0.0292
    # (r = 0.922), ka 0.411 (CV 64.1%); ka uncorrelated.
    etalcl + etalvc ~ c(0.0403, 0.0292, 0.0249)
    etalka ~ 0.411

    # Eq. 1 constant-CV error C = Cpred x (1 + eps); Table 4 sigma^2 = 0.105
    # (CV 32.4%), so the proportional SD is sqrt(0.105) = 0.324.
    propSd <- 0.324; label("Proportional residual error (SD, fraction)") # Table 4 final model: sigma^2 = 0.105 +/- 0.0060, CV 32.4%
  })
  model({
    # Linear-additive CL/F (mL/min) with exponential IIV on the
    # covariate-adjusted typical value; x 60 / 1000 converts mL/min to L/h.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 85) + e_wt_cl * (WT - 55)) *
      exp(etalcl) * 60 / 1000
    # Power-in-weight V/F plus the additive male shift (GEN = 1 - SEXF).
    vc <- (exp(lvc) * (WT / 55)^e_wt_vc + e_sex_vc * (1 - SEXF)) * exp(etalvc)
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and volume in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
