Paioni_2021_gentamicin <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in 109 hospitalised children aged 1 day to 14 years (median 29 days, median weight 4.23 kg) at University Children's Hospital Zurich receiving once-daily 5 or 7.5 mg/kg infusions over 2 or 30 min, fitted with the R package saemix (Paioni 2021). The authors parameterise the model by clearance, central volume V1, the distribution-phase rate constant lambda1 and an apparent peripheral volume V2' defined through the terminal rate constant lambda_z = CL / (V1 + V2'); these are converted here to the equivalent micro-constants. Clearance carries power effects of body weight (reference 4 kg), serum creatinine (reference 27 umol/L) and serum urea (reference 3 mmol/L); V1 and V2' carry power effects of body weight. Exponential (log-scale) residual error."
  reference <- paste(
    "Paioni P, Jaggi VF, Tilen R, Seiler M, Baumann P, Bram DS, Jetzer C,",
    "Haid RTU, Goetschi AN, Goers R, Muller D, Coman Schmid D,",
    "Meyer zu Schwabedissen HE, Rinn B, Berger C, Kramer SD.",
    "Gentamicin Population Pharmacokinetics in Pediatric Patients-A",
    "Prospective Study with Data Analysis Using the saemix Package in R.",
    "Pharmaceutics. 2021;13(10):1596.",
    "doi:10.3390/pharmaceutics13101596.",
    sep = " "
  )
  vignette <- "Paioni_2021_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gentamicin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (exponent 1.219), V1 (0.688) and V2' (0.974),",
        "all centred on 4 kg (Table 3 'Reference Value for Intercept'",
        "ln(4 kg)), close to the cohort median of 4.23 kg (Table 1, range",
        "2.5-35 kg). The paper does not say whether this is birth or",
        "current weight; with a median age of 29 days it is taken as the",
        "weight at the time of treatment."
      ),
      source_name = "body weight"
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL, (CREAT / 27)^-0.964 (Table 3). The reference",
        "27 umol/L is the assay's lower reporting limit: 73 of the 107",
        "subjects with a value were reported as '< 27 umol/L' (Table 1,",
        "range <27-72 umol/L), so the paper's 'close to the median'",
        "reference coincides with the censoring limit. The paper does not",
        "state how the censored values were coded; entering them at 27",
        "(zero log-difference) is the reading consistent with a reference",
        "at the limit, and the Discussion notes that urea 'may have",
        "compensated for the lacking serum creatinine values of < 27 uM'.",
        "Supply CREAT = 27 for a value reported below 27 umol/L; the",
        "fitted exponent is not informed by values below that limit.",
        "Missing values were entered at the reference (Section 2.3)."
      ),
      source_name = "serum creatinine"
    ),
    BUN = list(
      description = "Serum urea concentration (mmol/L; molar urea equals molar urea nitrogen).",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL, (BUN / 3)^-0.168 (Table 3), centred on 3 mmol/L",
        "(Table 1 median 3.20 mmol/L, range <1.8-17.3; 8 of 103 subjects",
        "reported as '< 1.8 mM'). The paper reports serum urea in mmol/L;",
        "convert a BUN in mg/dL with BUN_mmol = BUN_mgdL / 2.8. Missing",
        "values were entered at the reference (Section 2.3); the coding of",
        "the eight censored values is not stated."
      ),
      source_name = "serum urea"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on CL; reduced -2LL by less than 3.84 and was not retained",
        "(Section 3.3)."
      )
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Correlated as strongly with the CL random effect as body weight",
        "(adjusted r^2 0.576 vs 0.586) but not tested as an alternative",
        "because it is a calculated quantity (Section 3.3)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 109L,
    n_studies = 1L,
    n_observations = 310L,
    age_range = "0.9 days to 14.5 years (median 29.1 days, mean 213 days)",
    age_median = "29.1 days",
    weight_range = "2.5-35 kg",
    weight_median = "4.23 kg",
    sex_female_pct = 38.5,
    race_ethnicity = "Not reported.",
    ga_range = "29-42 weeks (median 39.3; N = 94); 11 born at < 37 weeks",
    disease_state = paste(
      "Hospitalised paediatric patients receiving gentamicin for at least",
      "48 h for suspected bacterial infection (43), suspected superimposed",
      "bacterial infection (14), proven bacterial infection (51; mostly",
      "pyelonephritis, neonatal sepsis and central-line bloodstream",
      "infection) or surgical prophylaxis (1) (Table 2). Patients with",
      "cystic fibrosis or on haemodialysis, peritoneal dialysis or",
      "haemofiltration were excluded; the Discussion notes no patients",
      "with impaired kidney function were included."
    ),
    dose_range = paste(
      "Once-daily IV gentamicin per hospital practice: 5 mg/kg at age",
      "< 7 days and 7.5 mg/kg at age >= 7 days, infused over 2 min or",
      "30 min depending on the unit."
    ),
    regions = "Single centre (University Children's Hospital Zurich, Switzerland).",
    notes = paste(
      "Prospective observational study, 1 October 2017 - 30 April 2019.",
      "Samples drawn before an infusion (trough) and about 30 min and 4 h",
      "after the end of the second or third infusion; 310 concentrations,",
      "71 below the LOQ of 0.3 mg/L. Concentrations from the Beckman",
      "GEN/DxC600 assay were transformed to the Abbott Alinity C scale",
      "(C_Alinity = 0.832 * C_DxC600 - 0.117 mg/L). Fitted with saemix 2.4",
      "in R 3.5.1 using an in-house BLQ treatment (residual set to zero when",
      "both predicted and observed are below the LOQ). Demographics from",
      "Table 1; parameter estimates from Table 3 (primary, non-bracketed",
      "values; the bracketed LOQ/2 sensitivity-analysis values are not",
      "used)."
    )
  )

  ini({
    # ----- Structural parameters (Paioni 2021 Table 3, final model) -----
    # Table 3 reports each intercept on the natural-log scale in PER-DAY
    # units (CL in L/d, lambda1 in 1/d); model() divides by 24 to work in
    # hours. Typical values at WT = 4 kg, CREAT = 27 umol/L, BUN = 3 mmol/L:
    #   CL = exp(2.267) = 9.65 L/d = 6.70 mL/min (Table 4 population
    #   median 6.88 mL/min at the cohort median weight 4.23 kg)
    #   lambda1 = exp(3.644) = 38.2 1/d = 1.59 1/h (Table 4 lambda1 1.59 1/h)
    lcl <- 2.267; label("Clearance at WT 4 kg, CREAT 27 umol/L, BUN 3 mmol/L (log of L/day)") # Paioni 2021 Table 3 'q1, Intercept (ln(CL, l x d-1))' = 2.267
    lvc <- 0.204; label("Central volume V1 at WT 4 kg (log of L)") # Paioni 2021 Table 3 'q4, Intercept (ln(V1, l))' = 0.204
    lvparea <- -0.174; label("Apparent peripheral volume V2' at WT 4 kg, with lambda_z = CL / (V1 + V2') (log of L)") # Paioni 2021 Table 3 'q2, Intercept (ln(V2', l))' = -0.174
    lalpha <- 3.644; label("Distribution-phase (fast) disposition rate constant lambda1 (log of 1/day)") # Paioni 2021 Table 3 'q3, Intercept (ln(l1, d-1))' = 3.644

    # ----- Covariate effects (Paioni 2021 Table 3) -----
    # Each enters as (ln(cov) - ln(ref)) on the log parameter, i.e. a power
    # of (cov / ref).
    e_wt_cl <- 1.219; label("Power exponent of (WT / 4 kg) on CL (unitless)") # Paioni 2021 Table 3 'q5, D ln(body weight, kg)' = 1.219
    e_creat_cl <- -0.964; label("Power exponent of (CREAT / 27 umol/L) on CL (unitless)") # Paioni 2021 Table 3 'q6, D ln(serum creatinine, uM)' = -0.964
    e_bun_cl <- -0.168; label("Power exponent of (BUN / 3 mmol/L) on CL (unitless)") # Paioni 2021 Table 3 'q7, D ln(serum urea, mM)' = -0.168
    e_wt_vparea <- 0.974; label("Power exponent of (WT / 4 kg) on V2' (unitless)") # Paioni 2021 Table 3 'q8, D ln(body weight, kg)' = 0.974
    e_wt_vc <- 0.688; label("Power exponent of (WT / 4 kg) on V1 (unitless)") # Paioni 2021 Table 3 'q9, D ln(body weight, kg)' = 0.688

    # ----- Inter-individual variability (Paioni 2021 Table 3) -----
    # Table 3 heading 'Variance of the random effects (inter-individual
    # variance)'; diagonal (no covariances reported). No random effect on
    # lambda1 (Section 3.2 and 3.3: omitted to avoid overparameterisation).
    # The variance reading is confirmed by Table 4: its individual-level
    # extremes of V1 and CL per 70 kg lie 0.30 and 0.77 log units beyond the
    # population-level extremes, i.e. 1.5 and 2.3 SDs under this reading but
    # more than 7 SDs if the printed values were SDs. The paper's own
    # Section 3.7 simulated C(30 min') intervals are narrower and match the
    # printed values used as SDs (ERRATUM in the source: its simulations
    # disagree with its own Table 3 scale); Table 3 is encoded verbatim as
    # variances. See the vignette section 'Errata in the source'.
    etalcl ~ 0.107 # Paioni 2021 Table 3 variance of random effects ln(CL) = 0.107
    etalvparea ~ 0.291 # Paioni 2021 Table 3 variance of random effects ln(V2') = 0.291
    etalvc ~ 0.0391 # Paioni 2021 Table 3 variance of random effects ln(V1) = 0.0391

    # ----- Residual error (Paioni 2021 Table 3) -----
    # saemix 'exponential' error model (Section 2.3): log(y) = log(f) + a*eps,
    # where saemix reports a as a standard deviation on the log scale.
    expSd <- 0.102; label("Exponential (log-scale) residual error SD") # Paioni 2021 Table 3 'Residual error' = 0.102
  })

  model({
    # Individual parameters in the source's parameterisation (Table 3
    # equations), converted from per-day to per-hour.
    cl <- exp(lcl + etalcl + e_wt_cl * log(WT / 4) + e_creat_cl * log(CREAT / 27) + e_bun_cl * log(BUN / 3)) / 24
    vc <- exp(lvc + etalvc + e_wt_vc * log(WT / 4))
    vparea <- exp(lvparea + etalvparea + e_wt_vparea * log(WT / 4))
    alpha <- exp(lalpha) / 24

    # Terminal rate constant (Section 2.3: lambda_z = CL / (V1 + V2')) and
    # the redistribution rate constant k21 = lambda1 * lambda_z * V1 / CL
    # (Equation 3). The remaining micro-constants follow from the
    # two-compartment identities lambda1 + lambda_z = kel + k12 + k21 and
    # lambda1 * lambda_z = kel * k21. The resulting ODE system has the same
    # eigenvalues and the same central-compartment solution as the paper's
    # closed-form Equations (1)-(2).
    kel <- cl / vc
    lambdaz <- cl / (vc + vparea)
    k21 <- alpha * lambdaz / kel
    k12 <- alpha + lambdaz - kel - k21
    vp <- vc * k12 / k21
    q <- k12 * vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
