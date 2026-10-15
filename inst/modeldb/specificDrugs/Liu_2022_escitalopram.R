Liu_2022_escitalopram <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption (ka",
    "fixed at 0.6 1/h) for oral escitalopram in 106 Chinese psychiatric",
    "inpatients (adolescents to older adults) sampled sparsely around trough",
    "during routine therapeutic drug monitoring (Liu 2022). Apparent clearance",
    "decreases linearly with age about the cohort median of 45 years and is",
    "multiplied by 0.847 in CYP2C19 intermediate and 0.479 in poor metabolizers",
    "relative to extensive metabolizers; exponential IIV on CL/F and V/F and a",
    "proportional residual error."
  )
  reference <- paste(
    "Liu S, Xiao T, Huang S, Li X, Kong W, Yang Y, Zhang Z, Ni X, Lu H,",
    "Zhang M, Shang D, Wen Y. Population pharmacokinetics model for",
    "escitalopram in Chinese psychiatric patients: effect of CYP2C19 and age.",
    "Front Pharmacol. 2022;13:964758. doi:10.3389/fphar.2022.964758.",
    sep = " "
  )
  vignette <- "Liu_2022_escitalopram"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Liu 2022 Methods 'Determination of escitalopram concentrations': serum
  # samples, oral conventional tablets (Results 'Demographic information').
  compartmentData <- list(
    depot = list(analyte = "escitalopram", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "escitalopram", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on CL/F centred on the cohort median of 45 years",
        "(Liu 2022 Eq 3 with the Table 1 median). Clearance DECREASES with",
        "age (see the e_age_cl comment in ini()). Observed range 12-83 years.",
        "Because the multiplier is linear it reaches zero at",
        "45 + 1/0.0077 = 175 years, far outside the data, so it stays",
        "positive over any realistic age."
      ),
      source_name = "AGE"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer indicator (1 = *1/*2 or *1/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Liu 2022 groups genotype into GENE = 1 (*1/*1, EM), 2 (*1/*2 or",
        "*1/*3, IM) and 3 (*2/*2 or *2/*3, PM) (Methods 'Determination of",
        "CYP2C19 genotype' and Eqs 5-7). CYP2C19_IM = 1 when GENE = 2.",
        "The reference category (both indicators 0) is the EM *1/*1 group;",
        "no *17 carriers were genotyped (Discussion). Cohort: EM 47, IM 49,",
        "PM 10 (Table 1)."
      ),
      source_name = "GENE"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer indicator (1 = *2/*2 or *2/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "CYP2C19_PM = 1 when the paper's GENE = 3; see the CYP2C19_IM entry for the grouping.",
      source_name = "GENE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and V/F by stepwise forward selection (Liu 2022 Methods, Eq 4 with COV = 1 for female) and not retained (Results 'PopPK model for escitalopram')."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Eq 3, median 61 kg) and not retained; missing weights were imputed to the median."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened (Eq 3, median 165 cm) and not retained; missing heights were imputed to the median."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Eq 3, median 22.4 kg/m^2) and not retained."
    ),
    SMOKE = list(
      description = "Current smoker",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained; only 3 of 106 patients smoked (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a liver-function covariate and not retained (median 17 U/L, Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened and not retained (median 40.4, Table 1; units not printed, g/L by magnitude)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as a renal-function covariate and not retained (median 67, Table 1; units not printed, umol/L by magnitude)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened with the other CYP2C19 inhibitors and inducers (Eq 4, 1 = co-administered during sampling) and not retained; 6 patients took omeprazole (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 106,
    n_studies = 1,
    n_observations = 337,
    age_range = "12-83 years",
    age_median = "45 years",
    weight_range = "37-97 kg",
    weight_median = "61 kg",
    sex_female_pct = 44.34,
    race_ethnicity = c(Asian = 100),
    disease_state = "Psychiatric inpatients receiving escitalopram with therapeutic drug monitoring",
    dose_range = "5-30 mg/day oral conventional tablets (5, 10, 15 or 20 mg once daily; 5 or 10 mg twice daily); median 10 mg/day",
    regions = "China (Affiliated Brain Hospital of Guangzhou Medical University, 2018-2021)",
    cyp2c19_phenotype = c(EM = 47, IM = 49, PM = 10),
    sampling_design = "Retrospective sparse TDM sampling, mostly around trough, at steady and non-steady state; patients with a single concentration were excluded.",
    notes = "Demographics from Liu 2022 Table 1; CYP2C19 allele and genotype frequencies from Table 2. Serum escitalopram by HPLC-MS/MS, linear range 3-300 ng/mL."
  )

  ini({
    # Structural parameters (Liu 2022 Table 4, final model). The reference
    # individual is a 45-year-old CYP2C19 extensive metabolizer: Eq 5 gives
    # CL = TVCL * theta_EM, and Table 4 estimates only theta_IM and theta_PM,
    # so theta_EM = 1.
    lka <- fixed(log(0.6)); label("Absorption rate constant (1/h)") # Table 4 'Ka (h-1) 0.6, FIX'; taken from Chen 2013 (Methods 'PopPK model development')
    lcl <- log(16.3); label("Apparent clearance CL/F for a 45-year-old CYP2C19 EM (L/h)") # Table 4 'CL/F (L/h) 16.3'
    lvc <- log(815); label("Apparent volume of distribution V/F (L)") # Table 4 'V/F (L) 815'

    # Covariate effects on CL/F (Table 4; Eqs 3 and 5-7).
    # Age: Table 4 prints theta_Age = 0.0077 and Eq 3 prints
    # CL = TVCL * (1 + theta * (AGE - 45)), which with a positive theta would
    # make clearance RISE with age. The Results text ('a decrease in CL/F of
    # escitalopram with increased patient age'), the Abstract and the typical
    # simulations in Figures 5 (65 years) and 6 (16 years) all have clearance
    # falling with age, so the effect is applied as (1 - e_age_cl * (AGE - 45)).
    e_age_cl <- 0.0077; label("Fractional decrease in CL/F per year of age above 45 years (1/year)") # Table 4 'theta Age 0.0077'; Eq 3; sign from Results text and Figures 5-6
    e_cyp2c19_im_cl <- 0.847; label("CL/F multiplier for CYP2C19 intermediate metabolizers (unitless)") # Table 4 'theta IM 0.847'; Eq 6
    e_cyp2c19_pm_cl <- 0.479; label("CL/F multiplier for CYP2C19 poor metabolizers (unitless)") # Table 4 'theta PM 0.479'; Eq 7

    # IIV (Eq 1, exponential). Table 4 'Random effect' values are NONMEM
    # OMEGA variances: sqrt(0.0877) = 29.6% CV-equivalent on CL/F and
    # sqrt(0.235) = 48.5% on V/F.
    etalcl ~ 0.0877 # Table 4 'Random effect CL/F 0.0877'
    etalvc ~ 0.235 # Table 4 'Random effect V/F 0.235'

    # Residual error (Eq 2, Y = F * (1 + eps1) + eps2). The additive term was
    # fixed to 0 (Table 4 'Additive error 0, FIX'), leaving a proportional
    # error with SIGMA variance 0.0287 -> SD = sqrt(0.0287) = 0.169.
    propSd <- sqrt(0.0287); label("Proportional residual error (fraction)") # Table 4 'Proportional error 0.0287' (variance)
  })
  model({
    # 1. Covariate terms on CL/F
    age_cl <- 1 - e_age_cl * (AGE - 45)
    cyp2c19_cl <- e_cyp2c19_im_cl^CYP2C19_IM * e_cyp2c19_pm_cl^CYP2C19_PM

    # 2. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * age_cl * cyp2c19_cl
    vc <- exp(lvc + etalvc)

    # 3. Micro-constant
    kel <- cl / vc

    # 4. ODEs
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # 5. Observation: dose mg / volume L = mg/L = ug/mL; x1000 gives ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
