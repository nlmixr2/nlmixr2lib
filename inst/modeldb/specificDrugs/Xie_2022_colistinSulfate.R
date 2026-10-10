Xie_2022_colistinSulfate <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous colistin sulfate in",
    "critically ill Chinese adults (Xie 2022; n = 20 ICU patients, 98 plasma",
    "concentrations collected after at least 72 h of therapy). Linear",
    "elimination from the central compartment with a 1-2 h intravenous",
    "infusion input. Cockcroft-Gault creatinine clearance enters clearance as",
    "a power function centred on 57.5 mL/min (exponent 0.353) and alanine",
    "aminotransferase enters the peripheral volume as a power function",
    "centred on 37 U/L (exponent 0.635). Inter-individual variability on CL",
    "and the central volume only; proportional residual error. Colistin",
    "sulfate is administered as the active drug and must not be confused",
    "with colistimethate sodium (CMS), the inactive prodrug. Dose unit",
    "conversion: 10,000 IU = 0.44 mg, so 1 million units (MU) = 44 mg."
  )
  reference <- paste(
    "Xie Y-L, Jin X, Yan S-S, Wu C-F, Xiang B-X, Wang H, Liang W, Yang B-C,",
    "Xiao X-F, Li Z-L, Pei Q, Zuo X-C, Peng Y (2022).",
    "Population pharmacokinetics of intravenous colistin sulfate and dosage",
    "optimization in critically ill patients.",
    "Front Pharmacol 13:967412.",
    "doi:10.3389/fphar.2022.967412.",
    sep = " "
  )
  vignette <- "Xie_2022_colistinSulfate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation",
        "(Methods, citing Cockcroft and Gault 1976), reported as RAW mL/min",
        "and NOT normalised to 1.73 m^2 body surface area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL centred on 57.5 mL/min, per Xie 2022 Eq. 3",
        "CL (L/h) = 1.50 x (CrCL/57.5)^0.353 x exp(eta_CL); the exponent is",
        "also the 'dCLdCrCL' row of Table 2. The centring constant 57.5",
        "appears ONLY inside Eq. 3 and is NOT the Table 1 cohort median of",
        "48.8 mL/min (range 6.5-193.8); the paper does not say which statistic",
        "it is. The paper does not print the body-weight convention used",
        "inside the Cockcroft-Gault equation. Two of 20 subjects received",
        "continuous renal replacement therapy (Table 1), and the paper gives",
        "no separate handling for them, so the model carries no explicit",
        "information about dialysis clearance. The paper's Monte Carlo",
        "simulations span 10-120 mL/min."
      ),
      source_name = "CrCL"
    ),
    ALT = list(
      description = "Serum alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the peripheral volume V2 centred on 37 U/L, per",
        "Xie 2022 Eq. 2 V2 (L) = 50.5 x (ALT/37)^0.635; the exponent is also",
        "the 'dV2dALT' row of Table 2. The Table 1 cohort median is 37.5 U/L",
        "(range 7-495). The paper's dosing simulations do not state an ALT",
        "value; V2 does not affect the steady-state average concentration."
      ),
      source_name = "ALT"
    )
  )

  # Screened in the stepwise covariate model and not retained (Results,
  # 'Population pharmacokinetic modeling': 'Age, sex, weight, and other
  # clinical variables had no statistically significant relationship with PK
  # parameters'). No point estimate is reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained. Table 1 median 60.5 (range 18-92) years.",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "unitless",
      type = "categorical",
      reference_category = "male",
      notes = "Screened and not retained. Table 1: 12 male / 8 female.",
      source_name = "sex"
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened and not retained. Table 1 median 55 (range 45-65) kg; the",
        "Discussion attributes the null finding to this narrow range."
      ),
      source_name = "weight"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained. Table 1 median 30.8 (range 22.2-43.3) g/L.",
      source_name = "ALB"
    )
  )

  compartmentData <- list(
    # Methods 'Quantification of colistin sulfate concentrations': colistin A
    # and colistin B in plasma were quantified by LC-MS/MS and 'the
    # concentrations of colistin A and B were added together and referred to
    # as colistin'. The analyte is active colistin, not a CMS-derived product.
    central = list(
      analyte = "colistin sulfate (sum of colistin A and colistin B)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "colistin sulfate (sum of colistin A and colistin B)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    n_observations = 98,
    age_range = "18-92 years (median 60.5)",
    weight_range = "45-65 kg (median 55)",
    sex_female_pct = 40,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Critically ill adults in the ICU receiving intravenous colistin sulfate",
      "for at least 72 h, predominantly for carbapenem-resistant Gram-negative",
      "infection (A. baumannii 50%, K. pneumoniae 30%, P. aeruginosa 10%,",
      "none cultured 20%); pulmonary infection 80%. APACHE II median 21.5",
      "(8-34), SOFA median 8 (3-18). Underlying liver disease 65%, chronic",
      "kidney disease 40%, sepsis 40%, solid-organ transplantation 20%."
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL median 48.8 mL/min (range 6.5-193.8); 2 of 20",
      "subjects on continuous renal replacement therapy."
    ),
    hepatic_function = "ALT median 37.5 U/L (range 7-495).",
    dose_range = paste(
      "1.0-2.0 MU daily (44-88 mg/day at 1 MU = 44 mg), most commonly",
      "1.5 MU/day, divided into 2 or 3 intravenous infusions of 1-2 h. Four",
      "subjects received a loading dose."
    ),
    co_medication = paste(
      "Concomitant meropenem 80%, ceftazidime/avibactam 15%,",
      "piperacillin/tazobactam, cefoperazone/sulbactam, linezolid and",
      "teicoplanin 10% each."
    ),
    regions = "China (single centre; Third Xiangya Hospital, Changsha)",
    notes = paste(
      "Single-centre study, May 2021 - April 2022. Demographics in Table 1.",
      "Plasma samples (38 troughs, 27 peaks, 33 random) were all collected",
      "after at least 72 h of therapy. Urinary recovery (6 subjects, 41 urine",
      "samples) was analysed by non-compartmental methods only and is not",
      "part of the model: median 10.05% of the dose (range 2.24-32.09%),",
      "median renal clearance 0.209 L/h. Modelled in Phoenix NLME 8.3.4 with",
      "FOCE-ELS."
    )
  )

  ini({
    # Table 2's 'CV (%)' column is the precision of each estimate (legend
    # 'percent confidence of variation'), quoted below as RSE.
    # Structural parameters: Xie 2022 Table 2 'Final model / Estimate',
    # identical to the typical values printed in Eqs 1-3.
    lcl <- log(1.50)
    label("Clearance at CRCL = 57.5 mL/min (L/h)") # Table 2 tvCL = 1.50 L/h (RSE 11.7%; Eq. 3)
    lvc <- log(16.1)
    label("Central volume of distribution, V (L)") # Table 2 tvV = 16.1 L (RSE 8.70%; Eq. 1)
    lvp <- log(50.5)
    label("Peripheral volume of distribution at ALT = 37 U/L, V2 (L)") # Table 2 tvV2 = 50.5 L (RSE 22.0%; Eq. 2)
    lq <- log(1.71)
    label("Intercompartmental clearance, CL2 (L/h)") # Table 2 tvCL2 = 1.71 L/h (RSE 22.0%); no IIV per Results text

    # Covariate effects: power functions (Methods 'continuous covariates ...
    # were modeled by using the power function'; Eqs 2-3).
    e_crcl_cl <- 0.353
    label("Power exponent on (CRCL / 57.5 mL/min) for CL (unitless)") # Table 2 dCLdCrCL = 0.353 (RSE 29.8%); Eq. 3
    e_alt_vp <- 0.635
    label("Power exponent on (ALT / 37 U/L) for V2 (unitless)") # Table 2 dV2dALT = 0.635 (RSE 16.6%); Eq. 2

    # Inter-individual variability: Table 2 reports omega^2 (variances) on V
    # and CL; CL2 and V2 carry no IIV.
    etalcl ~ 0.197 # Table 2 omega2 CL = 0.197 (RSE 29.2%)
    etalvc ~ 0.0267 # Table 2 omega2 V = 0.0267 (RSE 66.4%)

    # Residual error: proportional (Results); Table 2 'stdev0' is Phoenix
    # NLME's standard deviation of the proportional epsilon.
    propSd <- 0.228
    label("Proportional residual error (SD, fraction)") # Table 2 stdev0 = 0.228 (RSE 15.4%)
  })

  model({
    # Individual parameters (Eqs 1-3)
    cl <- exp(lcl + etalcl) * (CRCL / 57.5)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp) * (ALT / 37)^e_alt_vp
    q <- exp(lq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion dosed directly into central
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
