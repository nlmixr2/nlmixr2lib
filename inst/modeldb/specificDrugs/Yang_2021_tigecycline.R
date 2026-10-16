Yang_2021_tigecycline <- function() {
  description <- "One-compartment IV population PK model for tigecycline in Chinese adults with hospital-acquired pneumonia, complicated skin and skin-structure infection or complicated intra-abdominal infection, with blood urea nitrogen as a power covariate on clearance (Yang 2021)"
  reference <- "Yang T, Mei H, Wang J, Cai Y. Therapeutic Drug Monitoring of Tigecycline in 67 Infected Patients and a Population Pharmacokinetics/Microbiological Evaluation of A. baumannii Study. Front Microbiol. 2021;12:678165. doi:10.3389/fmicb.2021.678165"
  vignette <- "Yang_2021_tigecycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "SI units (mmol/L). Normalised to 9 mmol/L in the paper's final-model equation 'CL = 25.2 * (BUN/9)^-0.271' (Results, 'Construction of the Base Model and Full Regression Model'); 9 mmol/L is the upper bound of the paper's compensatory renal-insufficiency grade (Table 5 footnote), not a stated cohort median. Enters CL as a power term with a negative exponent, so CL falls as BUN rises. Observed range 1.8-48.3 mmol/L (Table 5 strata). Baseline vs time-varying status not stated.",
      source_name = "BUN"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Significant on CL in the univariable and multivariable steps but removed at backward elimination (Supplement, 'Results of the full regression model (FRM) establishment'); no estimate reported."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a candidate covariate (Methods, 'Modeling of the PPK/PD of TGC') and not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a candidate covariate (Methods) and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as a candidate covariate (Methods, 'Cr') and not retained; units not printed."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened as a candidate covariate (Methods) and not retained; units not printed."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as a candidate covariate (Methods, 'TBIL') and not retained; units not printed. Direct bilirubin, ALT, AST, GGT and ALP were screened the same way."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 67L,
    n_studies = 1L,
    n_observations = 186L,
    age_mean = "65.5 years (SD 23)",
    weight_mean = "64.83 kg (SD 14.96)",
    sex_female_pct = 25.4,
    race_ethnicity = "Chinese (single-centre study in Beijing)",
    disease_state = "Adults (> 18 years) with hospital-acquired pneumonia (HAP), complicated skin and skin-structure infection (cSSSI) or complicated intra-abdominal infection (cIAI) treated with tigecycline; infection site in the 58 patients with outcome data: lung 55%, lung plus another site 17%, abdomen 7%, other 21%. Acinetobacter baumannii was isolated in 41% of those patients (alone or with K. pneumoniae / E. coli).",
    dose_range = "IV tigecycline at steady state: 100 mg loading dose then 50 mg q12h (n = 24), 100 mg q12h (n = 29), or 50 mg q12h (n = 14)",
    regions = "China (PLA General Hospital, Beijing)",
    renal_function = "Blood urea nitrogen 1.8-48.3 mmol/L across the simulation strata of Table 5 (normal 1.8-7.5, compensatory renal insufficiency up to 9, decompensated > 9, renal failure > 20 mmol/L)",
    notes = "Demographics from Yang 2021 Results ('Profiles of Study Patients') and Table 1 (50 male / 17 female; per-regimen counts 16/8, 22/7 and 12/2 male/female). Prospective therapeutic-drug-monitoring study, October 2018 - March 2019. Samples were taken at steady state at one or more of 0, 0.5, 1.5, 3, 6, 9 and 12 h after a dose, mostly troughs; 186 samples in total, one below the LLOQ of 50 ng/mL. Model fitted in NONMEM 6; bootstrap with 1000 replicates converged in 100%."
  )

  ini({
    # Structural parameters: Yang 2021 Supplementary Table S10, 'Population
    # pharmacokinetic parameters of TGC final regression model'. The same
    # values appear in the main-text equation 'CL = 25.2 * (BUN/9)^-0.271'.
    lcl <- log(25.2); label("Clearance at BUN = 9 mmol/L (L/h)") # Table S10: CL, L/h = 25.2 (RSE 7.9%)
    lvc <- log(417); label("Volume of distribution (L)") # Table S10: V, L = 417 (RSE 33.3%)

    # Covariate effect: Table S10 'theta BUN (CL)' = -0.271, power exponent on
    # (BUN / 9).
    e_bun_cl <- -0.271; label("Power exponent on (BUN / 9 mmol/L) for CL (unitless)") # Table S10: theta BUN (CL) = -0.271 (RSE 43.2%)

    # Inter-individual variability on CL only (Table S10 'Inter-trial
    # variability'). The table prints 0.587 without saying whether it is
    # omega (SD) or omega^2 (variance); it is read as the SD, so the
    # variance is 0.587^2. See the vignette for the variance reading.
    etalcl ~ 0.344569 # 0.587^2; Table S10: eta (CL) = 0.587 (RSE 11.3%), read as an SD

    # Residual error: proportional (Supplement, 'Base model': 'The proportional
    # without ETA residual-error model was taken as the final base model').
    # Table S10 prints 0.576; read as the SD, consistent with the eta reading.
    propSd <- 0.576; label("Proportional residual error (fraction)") # Table S10: epsilon (prop) = 0.576 (RSE 10.5%), read as an SD
  })
  model({
    # Individual parameters
    cl <- exp(lcl + etalcl) * (BUN / 9)^e_bun_cl
    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment model; tigecycline is given by IV infusion into central
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> Cc in mg/L (= ug/mL, the paper's unit)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
