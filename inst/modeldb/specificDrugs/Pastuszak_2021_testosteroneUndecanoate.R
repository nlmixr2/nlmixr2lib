Pastuszak_2021_testosteroneUndecanoate <- function() {
  description <- "One-compartment population PK model for total serum testosterone after intramuscular testosterone undecanoate 750 mg in adult men with hypogonadism, with first-order absorption from the injection depot, linear apparent clearance, and an apparent zero-order endogenous testosterone input into the central compartment that is acutely suppressed from its predose rate after the first injection and then recovers gradually to a steady-state rate. Body weight scales CL/F (fixed allometric 0.75), V/F (fixed 1) and Ka (estimated -1.83); baseline sex hormone-binding globulin scales CL/F."
  reference <- paste(
    "Pastuszak AW, Bush M, Curd L, Vijayan S, Priestley T, Xiang Q, Hu Y.",
    "Population Pharmacokinetic Modeling and Simulations to Evaluate a Potential Dose Regimen of",
    "Testosterone Undecanoate in Hypogonadal Males.",
    "J Clin Pharmacol 2021;61(12):1618-1625.",
    "doi:10.1002/jcph.1939.",
    sep = " "
  )
  vignette <- "Pastuszak_2021_testosteroneUndecanoate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/dL")

  # Amounts are mg and the volume is L, so central / vc is mg/L; 1 mg/L =
  # 1e5 ng/dL, the unit total testosterone is reported in throughout the paper.
  compartmentData <- list(
    depot = list(analyte = "testosterone undecanoate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "testosterone (total)", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value, time-fixed per subject (Methods: 'Baseline demographics of age, albumin level, body weight, ... were evaluated'). Power scaling referenced to 101 kg (Table 2 footnote b) on CL/F (exponent 0.75 fixed), V/F (1 fixed) and Ka (-1.83 estimated). Table 1 PK population mean 101.6 (SD 17.4) kg; the original study PK population required weight >= 65 kg.",
      source_name = "Weight"
    ),
    SEX_HORMONE_BINDING_GLOBULIN = list(
      description = "Baseline serum sex hormone-binding globulin concentration.",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value, time-fixed per subject. Power scaling referenced to 20 nmol/L (Table 2 footnote b) on CL/F, exponent -0.219: lower SHBG, higher CL/F. Table 1 PK population mean 20.8 (SD 8.5) nmol/L.",
      source_name = "SHBG"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and V/F (Methods, Population PK Modeling) but not retained in the final model."
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and V/F but not retained. Table 1 mean 4.2 (SD 0.3) g/dL = 42 g/L."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F, V/F and the rate and extent of absorption but not retained; body weight was preferred (Discussion). Table 1 mean 32.0 (SD 5.2) kg/m^2."
    ),
    RACE_BLACK = list(
      description = "Black race indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "non-Black",
      notes = "Race (White / Black / Other) was screened on CL/F and V/F but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 130L,
    n_studies = 1L,
    n_observations = 2360L,
    age_range = "24-75 years (enrolled); 30-75 years (PK population)",
    age_mean = "54.2 years (enrolled); 54.6 years (PK population)",
    weight_mean = "101.2 kg (SD 18.0) enrolled; 101.6 kg (SD 17.4) PK population",
    sex_female_pct = 0,
    race_ethnicity = c(White = 74.6, Black = 12.3, Other = 13.1),
    disease_state = "adult men with hypogonadism (primary or hypogonadotropic)",
    dose_range = "testosterone undecanoate 750 mg in 3 mL castor oil / benzyl benzoate intramuscularly into the buttock at baseline, week 4 and then every 10 weeks for up to 84 weeks",
    regions = "United States (31 sites)",
    shbg_mean = "20.9 nmol/L (SD 8.9) enrolled; 20.8 nmol/L (SD 8.5) PK population",
    notes = "Phase 3 single-arm open-label Study IP157-001 Part C (NCT00467870); Table 1. Dense sampling across injection intervals 3 and 4 (weeks 14-34) with sparser sampling elsewhere; testosterone measured by LC-MS/MS (LLOQ 20.0 ng/dL). Data after a crossover to subcutaneous dosing (Part D) were excluded. 117 of the 130 men formed the original study PK population used as the simulation template."
  )

  ini({
    lka <- log(0.001); label("First-order absorption rate constant from the intramuscular depot at 101 kg (Ka, 1/h)") # Table 2, Ka 0.001 (units printed 'L/h')
    lcl <- log(197); label("Apparent clearance of total testosterone at 101 kg and SHBG 20 nmol/L (CL/F, L/h)") # Table 2, CL/F 197 L/h
    lvc <- log(12300); label("Apparent volume of distribution at 101 kg (V/F, L)") # Table 2, V/F 12 300 L

    lkin_te <- log(0.445); label("Apparent baseline (predose) endogenous testosterone production rate (Rb, mg/h)") # Table 2, Rb 0.445 mg/h
    lkin_te_ss <- log(0.572); label("Apparent steady-state endogenous testosterone production rate (Rss, mg/h)") # Table 2, Rss 0.572 mg/h
    lk_suppression_te <- fixed(log(10)); label("First-order rate constant for acute suppression of endogenous testosterone production (K1, 1/h)") # Table 2, K1 10 FIXED (units printed 'L/h')
    lk_recovery_te <- log(0.000727); label("First-order rate constant for gradual recovery of endogenous testosterone production (K2, 1/h)") # Table 2, K2 0.000727 (units printed 'L/h')

    e_wt_cl <- fixed(0.75); label("Power exponent of (WT / 101 kg) on CL/F (unitless)") # Table 2, 'Weight, CL/F (power)' 0.75 FIXED
    e_sex_hormone_binding_globulin_cl <- -0.219; label("Power exponent of (SEX_HORMONE_BINDING_GLOBULIN / 20 nmol/L) on CL/F (unitless)") # Table 2, 'SHBG, CL/F (power)' -0.219
    e_wt_ka <- -1.83; label("Power exponent of (WT / 101 kg) on Ka (unitless)") # Table 2, 'Weight, Ka (power)' -1.83
    e_wt_vc <- fixed(1); label("Power exponent of (WT / 101 kg) on V/F (unitless)") # Table 2, 'Weight, V/F (power)' 1.0 FIXED

    etalcl ~ 0.0392207 # Table 2, IIV CL/F 20.0 %CV; log(0.200^2 + 1)
    etalka ~ 0.200361 # Table 2, IIV Ka 47.1 %CV; log(0.471^2 + 1)

    propSd <- 0.188; label("Proportional residual error (fraction)") # Table 2, proportional residual error 18.8%
    addSd <- 56.8; label("Additive residual error (ng/dL)") # Table 2, additive residual error SE 56.8 ng/dL
  })

  model({
    ka <- exp(lka + etalka) * (WT / 101)^e_wt_ka
    cl <- exp(lcl + etalcl) * (WT / 101)^e_wt_cl * (SEX_HORMONE_BINDING_GLOBULIN / 20)^e_sex_hormone_binding_globulin_cl
    vc <- exp(lvc) * (WT / 101)^e_wt_vc
    kel <- cl / vc

    kin_te <- exp(lkin_te)
    kin_te_ss <- exp(lkin_te_ss)
    k_suppression_te <- exp(lk_suppression_te)
    k_recovery_te <- exp(lk_recovery_te)

    # Apparent endogenous testosterone input (Results, equation under Figure 1):
    # R = Rb * exp(-K1 * TIME) + Rss * (1 - exp(-K2 * TIME)), with TIME the time
    # since the first injection.
    rate_endogenous_te <- kin_te * exp(-k_suppression_te * t) + kin_te_ss * (1 - exp(-k_recovery_te * t))

    # Predose endogenous steady state, Rb / CL in concentration terms.
    central(0) <- kin_te / kel

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central + rate_endogenous_te

    Cc <- central / vc * 1e5
    Cc ~ prop(propSd) + add(addSd)
  })
}
