Nakao_2019_afatinib <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral afatinib in Japanese adults with EGFR mutation-positive non-small cell lung cancer, with centred-linear AST and creatinine-clearance effects on CL/F and centred-linear BMI and age effects on V/F (Nakao 2019)"
  reference <- "Nakao K, Kobuchi S, Marutani S, Iwazaki A, Tamiya A, Isa S, Okishio K, Kanazu M, Tamiya M, Hirashima T, Imai K, Sakaeda T, Atagi S. Population pharmacokinetics of afatinib and exposure-safety relationships in Japanese patients with EGFR mutation-positive non-small cell lung cancer. Sci Rep. 2019;9:18202. doi:10.1038/s41598-019-54804-9"
  vignette <- "Nakao_2019_afatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) in the source. Enters CL/F as the centred-linear factor (1 + (AST - 25.3) * e_ast_cl) (Nakao 2019 Table 2 CL/F equation). Cohort mean 25.6 IU/L, range 13-65 (Table 1). The factor reaches zero at AST ~ 87.8 IU/L, above the observed range; do not extrapolate.",
      source_name = "AST"
    ),
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalized (Nakao 2019 Table 1 footnote and Methods 'Individual Ccr values were determined using the Cockcroft-Gault equation'). Enters CL/F as (1 + (CRCL - 79.9) * e_crcl_cl) (Table 2). Cohort mean 80.8 mL/min, range 42.3-131.8 (Table 1).",
      source_name = "Ccr"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters V/F as (1 + (BMI - 21.8) * e_bmi_vc) (Nakao 2019 Table 2 V/F equation). Cohort mean 21.9 kg/m^2, range 15.2-28.1 (Table 1).",
      source_name = "BMI"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters V/F as (1 + (AGE - 66.7) * e_age_vc) (Nakao 2019 Table 2 V/F equation). Cohort mean 66.8 years, range 45-86 (Table 1).",
      source_name = "Age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained (Nakao 2019 Methods 'Development of population pharmacokinetics model').",
      source_name = "Weight"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained (Nakao 2019 Methods).",
      source_name = "Height"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained (Nakao 2019 Methods; Discussion notes AST but not ALT was retained).",
      source_name = "ALT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained (Nakao 2019 Methods).",
      source_name = "Cre"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "afatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "afatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_observations = 354L,
    age_range = "45-86 years",
    age_mean = "66.8 years",
    weight_range = "35.5-79.1 kg",
    weight_mean = "53.8 kg",
    bmi_range = "15.2-28.1 kg/m^2",
    sex_female_pct = 67.6,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Advanced non-small cell lung cancer (all adenocarcinoma) harbouring an EGFR-activating mutation (Ex19 del 18/34, L858R 9/34, uncommon/compound mutations 7/34); 24/34 had prior chemotherapy and 21/34 a prior EGFR TKI (gefitinib or erlotinib).",
    dose_range = "Afatinib 40 mg orally once daily (dose escalation to 50 mg permitted at physician discretion; administered fasted).",
    regions = "Japan (Kinki-Chuo Chest Medical Center and Osaka Habikino Medical Center, August 2014 - May 2016; UMIN000014181).",
    renal_function = "Cockcroft-Gault CrCl mean 80.8 mL/min, range 42.3-131.8; 10 patients below the normal limit.",
    hepatic_function = "AST mean 25.6 IU/L (range 13-65), ALT mean 19.7 IU/L (range 5-77); 4 patients with AST/ALT above the normal limit.",
    notes = "Nakao 2019 Table 1 (baseline demographics; 11 male / 23 female, ECOG PS 0-1 in 28/34). Sampling at 0.5-1, 2-3, 4-6, 8-12 and 24 h on day 1 and 0-24 h on day 8; afatinib assayed by HPLC. Phoenix NLME 7.0, FOCE-ELS."
  )

  ini({
    lka <- log(0.60); label("Absorption rate constant ka (1/h)") # Nakao 2019 Table 2, 'ka (1/h)' = 0.60
    lcl <- log(20.0); label("Apparent clearance CL/F for the typical patient (L/h)") # Nakao 2019 Table 2, 'theta CL (L/h)' = 20.0
    lvc <- log(795.8); label("Apparent volume of distribution V/F for the typical patient (L)") # Nakao 2019 Table 2, 'theta V (L)' = 795.8

    e_crcl_cl <- 0.0013; label("Centred-linear slope of creatinine clearance on CL/F (min/mL)") # Nakao 2019 Table 2, 'theta Ccr (min/mL)' = 0.0013
    e_ast_cl <- -0.016; label("Centred-linear slope of AST on CL/F (L/IU)") # Nakao 2019 Table 2, 'theta AST (L/IU)' = -0.016
    e_bmi_vc <- 0.019; label("Centred-linear slope of BMI on V/F (m^2/kg)") # Nakao 2019 Table 2, 'theta BMI (m2/kg)' = 0.019
    e_age_vc <- -0.004; label("Centred-linear slope of age on V/F (1/year)") # Nakao 2019 Table 2, 'theta Age' = -0.004

    # IIV: Table 2 reports omega as a percentage; interpreted as a log-normal
    # CV%, omega^2 = log(1 + CV^2) (see vignette Assumptions for the evidence).
    etalka ~ 0.62221 # Nakao 2019 Table 2, 'omega ka (%)' = 92.9 -> log(1 + 0.929^2)
    etalcl ~ 0.46169 # Nakao 2019 Table 2, 'omega CL/F (%)' = 76.6 -> log(1 + 0.766^2)
    etalvc ~ 0.24591 # Nakao 2019 Table 2, 'omega V/F (%)' = 52.8 -> log(1 + 0.528^2)

    propSd <- 0.317; label("Proportional residual error (fraction)") # Nakao 2019 Table 2, 'sigma (%)' = 31.7
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (1 + (CRCL - 79.9) * e_crcl_cl) * (1 + (AST - 25.3) * e_ast_cl)
    vc <- exp(lvc + etalvc) * (1 + (BMI - 21.8) * e_bmi_vc) * (1 + (AGE - 66.7) * e_age_vc)

    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L = ug/mL; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
