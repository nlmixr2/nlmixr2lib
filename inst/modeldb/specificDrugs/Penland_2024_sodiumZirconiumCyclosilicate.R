Penland_2024_sodiumZirconiumCyclosilicate <- function() {
  description <- paste(
    "Population pharmacodynamic (K-PD) indirect-response model for the",
    "effect of oral sodium zirconium cyclosilicate (SZC), a non-absorbed",
    "potassium binder, on serum potassium in adults with hyperkalemia",
    "(Penland 2024; pooled phase II/III studies). SZC doses enter a virtual",
    "exposure compartment with first-order elimination (mean transit time",
    "137 h); the exposure is scaled so that a steady-state 1 g once-daily",
    "regimen has a mean 24-h exposure of 1, which makes EC50 a once-daily",
    "dose (g). Serum potassium follows an indirect-response model whose",
    "production rate kin = baseline * kout is multiplied by a constant",
    "on-treatment placebo factor and inhibited by a sigmoid Imax function of",
    "the SZC exposure. Final full-covariate model: baseline serum potassium,",
    "age, body weight, eGFR and race on EC50; baseline serum potassium, eGFR,",
    "heart failure and diabetes on the placebo factor; sex on kout.",
    "Log-normal IIV on the placebo factor and EC50; additive residual error."
  )
  reference <- paste(
    "Penland RC, Astrand M, Boulton DW, Nagard M. Population",
    "Pharmacodynamic Dose-Response Analysis of Serum Potassium Following",
    "Dosing with Sodium Zirconium Cyclosilicate. Clin Pharmacokinet.",
    "2024;63(4):551-560. doi:10.1007/s40262-024-01360-9.",
    "Model equations and covariate functional forms are in the",
    "Supplementary Methods (Electronic Supplementary Material).",
    sep = " "
  )
  vignette <- "Penland_2024_sodiumZirconiumCyclosilicate"
  units <- list(
    time = "h",
    dosing = "g",
    concentration = "mmol/L (serum potassium)"
  )

  compartmentData <- list(
    depot_kpd = list(
      analyte = "sodium zirconium cyclosilicate",
      units = "g",
      specimen = "not applicable",
      verified = TRUE
    ),
    serum_k = list(analyte = "potassium", units = "mmol/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    POT = list(
      description = "Baseline serum potassium at the time of the first SZC dose",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline. Serves both as the initial condition",
        "serum_k(0) = POT (Supplementary Methods: SK(0) = BL at t = 0, the",
        "time of first dose) and as a power covariate (POT / 5.4) on EC50",
        "(exponent -1.04) and on the placebo factor (exponent 0.44). The",
        "centring value 5.4 mmol/L is the paper's reference patient."
      ),
      source_name = "BL"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate (AGE / 65)^-0.486 on EC50; 65 years is the paper's reference patient.",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate (WT / 85)^0.515 on EC50; 85 kg is the paper's reference patient.",
      source_name = "Body weight"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the MDRD equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "MDRD eGFR (Table 1 footnote a). Power covariate (CRCL / 45) on EC50",
        "(exponent 0.316) and on the placebo factor (exponent 0.0371);",
        "45 mL/min/1.73 m^2 is the paper's reference patient."
      ),
      source_name = "eGFR"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White)",
      notes = "Proportional shift exp(-0.325 * RACE_BLACK) on EC50; White race is the reference.",
      source_name = "Black/African American"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White)",
      notes = "Proportional shift exp(-0.134 * RACE_ASIAN) on EC50; White race is the reference.",
      source_name = "Asian"
    ),
    RACE_OTHER = list(
      description = "Other race indicator (not White, Black/African American or Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White)",
      notes = "Proportional shift exp(-0.122 * RACE_OTHER) on EC50; White race is the reference (Table 1 footnote b defines 'Other').",
      source_name = "Other race"
    ),
    DIS_CHF = list(
      description = "Congestive heart failure (medical history) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no congestive heart failure)",
      notes = "Proportional shift exp(0.00121 * DIS_CHF) on the placebo factor; 1 = medical history of congestive heart failure (Table 1).",
      source_name = "Heart failure status"
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus (medical history) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes)",
      notes = "Proportional shift exp(-0.0133 * DIS_DIAB) on the placebo factor.",
      source_name = "Diabetes status"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Proportional shift exp(0.254 * SEXF) on kout. The paper does not",
        "print the coding of its sex covariate; it is taken as 1 = female",
        "because the paper's reference patient is male and the categorical",
        "form gives the typical value at covariate = 0. kout cancels from",
        "the steady state, so the coding only affects the speed of the",
        "serum potassium response."
      ),
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2369L,
    n_studies = 7L,
    n_observations = 25764L,
    age_range = "mean (SD) 66.0 (12) years",
    weight_range = "mean (SD) 87.6 (22) kg (Table 1 total column)",
    sex_female_pct = 39,
    race_ethnicity = c(White = 68, Black = 9, Asian = 21, Other = 2),
    disease_state = "hyperkalemia (screening serum potassium > 5.0 mmol/L); 34% congestive heart failure, 62% diabetes",
    potassium_range = "baseline mean (SD) 5.49 (0.43) mmol/L",
    renal_function = "MDRD eGFR mean (SD) 45.9 (29), range 3.46-190 mL/min/1.73 m^2",
    dose_range = "correction phase 0-15 g three times daily for 48-72 h; maintenance phase 0-15 g once daily up to 28 days",
    regions = "global, including three Asian studies (D9480C00002, D9482C00001, D9482C00002)",
    notes = "Studies ZS-002, ZS-003, ZS-004, ZS-005, D9480C00002, D9482C00001 and D9482C00002 (Table 1; Supplementary Table S1)."
  )

  ini({
    # Virtual SZC exposure (K-PD) compartment
    lkel <- log(1 / 137); label("First-order elimination rate of the virtual SZC exposure compartment, 1/MTT (1/h)") # Table 2 full model 'SZC MTT, h' = 137; kel = 1/MTT

    # Serum potassium indirect-response model
    lkout <- log(0.548); label("First-order serum potassium loss rate kout (1/h)") # Table 2 full model 'Serum K+ Kout, 1/h' = 0.548
    lplac <- log(0.923); label("Constant on-treatment placebo factor multiplying kin (unitless)") # Table 2 full model 'Placebo effect on Kin' = 0.923
    limax <- log(0.633); label("Maximum fractional inhibition of kin by SZC, Imax (fraction)") # Table 2 full model 'Emax, % inhibition' = 63.3
    lec50 <- log(32.8); label("SZC exposure at half-maximal inhibition, EC50, expressed as a once-daily dose (g)") # Table 2 full model 'EC50, g' = 32.8
    lhill <- log(1.36); label("Hill coefficient of the SZC inhibition (unitless)") # Table 2 full model 'Hill coefficient' = 1.36

    # Covariate effects on EC50 (power for continuous, exp() for categorical)
    e_pot_ec50 <- -1.04; label("Power exponent of (POT / 5.4) on EC50 (unitless)") # Table 2 'BL covariate on EC50' = -1.04
    e_age_ec50 <- -0.486; label("Power exponent of (AGE / 65) on EC50 (unitless)") # Table 2 'Age covariate on EC50' = -0.486
    e_wt_ec50 <- 0.515; label("Power exponent of (WT / 85) on EC50 (unitless)") # Table 2 'Body weight covariate on EC50' = 0.515
    e_crcl_ec50 <- 0.316; label("Power exponent of (CRCL / 45) on EC50 (unitless)") # Table 2 'eGFR covariate on EC50' = 0.316
    e_race_black_ec50 <- -0.325; label("Log proportional shift in EC50 for Black/African American race (unitless)") # Table 2 'Black/African American covariate on EC50' = -0.325
    e_race_asian_ec50 <- -0.134; label("Log proportional shift in EC50 for Asian race (unitless)") # Table 2 'Asian covariate on EC50' = -0.134
    e_race_other_ec50 <- -0.122; label("Log proportional shift in EC50 for Other race (unitless)") # Table 2 'Other race covariate on EC50' = -0.122

    # Covariate effects on the placebo factor
    e_pot_plac <- 0.44; label("Power exponent of (POT / 5.4) on the placebo factor (unitless)") # Table 2 'BL covariate on placebo' = 0.44
    e_crcl_plac <- 0.0371; label("Power exponent of (CRCL / 45) on the placebo factor (unitless)") # Table 2 'eGFR covariate on placebo' = 0.0371
    e_dis_chf_plac <- 0.00121; label("Log proportional shift in the placebo factor for congestive heart failure (unitless)") # Table 2 'Heart failure status covariate on placebo' = 0.00121
    e_dis_diab_plac <- -0.0133; label("Log proportional shift in the placebo factor for diabetes (unitless)") # Table 2 'Diabetes status covariate on placebo' = -0.0133

    # Covariate effect on kout
    e_sexf_kout <- 0.254; label("Log proportional shift in kout for female sex (unitless)") # Table 2 'Sex covariate on Kout' = 0.254

    # Between-patient variability (log-normal; Supplementary Methods)
    etalplac ~ 0.00296 # Table 2 full model 'BPV of placebo variance' = 0.00296
    etalec50 ~ 0.171 # Table 2 full model 'BPV of EC50 variance' = 0.171

    # Additive residual error on serum potassium (Supplementary Methods)
    addSd <- sqrt(0.118); label("Additive residual error on serum potassium (mmol/L)") # Table 2 full model 'Residual error variance' = 0.118; SD = sqrt(0.118)
  })
  model({
    kel <- exp(lkel)

    ec50 <- exp(lec50 + etalec50) *
      (POT / 5.4)^e_pot_ec50 *
      (AGE / 65)^e_age_ec50 *
      (WT / 85)^e_wt_ec50 *
      (CRCL / 45)^e_crcl_ec50 *
      exp(e_race_black_ec50 * RACE_BLACK + e_race_asian_ec50 * RACE_ASIAN + e_race_other_ec50 * RACE_OTHER)

    plac <- exp(lplac + etalplac) *
      (POT / 5.4)^e_pot_plac *
      (CRCL / 45)^e_crcl_plac *
      exp(e_dis_chf_plac * DIS_CHF + e_dis_diab_plac * DIS_DIAB)

    kout <- exp(lkout) * exp(e_sexf_kout * SEXF)
    imax <- exp(limax)
    hill <- exp(lhill)

    # Virtual SZC exposure: depot_kpd holds grams of SZC. At steady state on
    # D g once daily the mean amount is D / (kel * 24 h), so multiplying by
    # kel * 24 scales the exposure to the once-daily dose equivalent (g);
    # 1 g QD gives a mean 24-h exposure of 1, as in the paper.
    d/dt(depot_kpd) <- -kel * depot_kpd
    exposure <- depot_kpd * kel * 24

    inh <- imax * exposure^hill / (exposure^hill + ec50^hill)

    # Supplementary Methods: Kin = SK(0) * Kout * placebo * (1 - inhibition),
    # with the placebo factor applied from t = 0, the time of the first dose
    kin <- POT * kout * plac * (1 - inh)
    d/dt(serum_k) <- kin - kout * serum_k
    serum_k(0) <- POT

    serum_k ~ add(addSd)
  })
}
