Ge_2020_metoclopramide <- function() {
  description <- "Two-compartment population PK model with first-order absorption for metoclopramide in infants, children, and adolescents (Ge 2020), fit simultaneously to opportunistic plasma concentrations after intravenous (bolus or infusion) and enteral (oral, nasogastric, nasojejunal, gastrostomy) dosing. Body weight is the only covariate: clearance and intercompartmental clearance scale allometrically on WT/70 with a fixed exponent of 0.75, and both volumes scale linearly (fixed exponent 1). A single first-order absorption rate and a single bioavailability apply to all enteral routes; interindividual variability is estimated on clearance only, with a proportional residual error."
  reference <- "Ge S, Mendley SR, Gerhart JG, Melloni C, Hornik CP, Sullivan JE, Atz A, Delmore P, Tremoulet A, Harper B, Payne E, Lin S, Erinjeri J, Cohen-Wolkowiez M, Gonzalez D; Best Pharmaceuticals for Children Act - Pediatric Trials Network Steering Committee. Population Pharmacokinetics of Metoclopramide in Infants, Children, and Adolescents. Clin Transl Sci. 2020;13(6):1189-1198. doi:10.1111/cts.12803"
  vignette <- "Ge_2020_metoclopramide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "metoclopramide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "metoclopramide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "metoclopramide", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: the most recent measured weight for each record (Ge 2020 Methods 'Covariate selection'). Median (5th-95th percentile) 23.5 (2.6-98.1) kg at the first recorded dose (Table 1). Included a priori: CL and Q scale as (WT/70)^0.75 and Vc and Vp as (WT/70)^1, both exponents fixed (Eqs. 6-9).",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Screened as a continuous covariate, as PNA group A (<1 month, 1 month-<2 years, 2-<12 years, 12-<17 years, 17-<21 years) and as PNA group B (<1 month vs >=1 month), and in linear, power and sigmoidal Emax maturation functions on CL. PNA group B entered in the forward step (dOFV -4.3; CL about 30% lower below 1 month) but was removed in backward elimination (Ge 2020 Results 'PopPK model development and evaluation', Discussion). Median 8.9 years (5th-95th percentile 0.02-18.5; Table 1)."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened on CL but not retained (Ge 2020 Methods 'Covariate selection'). Collected only for the 17 infants <120 days PNA: median 37.9 (5th-95th percentile 29.5-40.5) weeks (Table 1)."
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      notes = "Screened on CL (linear, power and sigmoidal Emax maturation) but not retained (Ge 2020 Methods 'Covariate selection', Discussion)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "The only covariate significant on CL after backward elimination (dOFV -7.0; power exponent -0.27), but excluded from the final model: it did not reduce IIV in CL, destabilised the model, failed the covariance step, and lost significance after dropping one subject with SCr 4.3 mg/dL (Ge 2020 Results, Figure S1, Discussion). Available for 31 of 50 subjects; median 0.5 (0.2-1.1) (Table 1, printed as mg/mL, which is a typographical slip for mg/dL)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Schwartz equation)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 109.7 (5.0-225.1) mL/min/1.73 m^2 in 31 subjects (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 3.5 (2.4-4.4) g/dL in 17 subjects (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 39.5 (17.6-116.0) U/L in 14 subjects (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 39.5 (17.1-154.1) U/L in 12 subjects (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 0.7 (0.2-15.2) mg/dL in 15 subjects (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'). Median 19.0 (9.7-35.8) (Table 1)."
    ),
    DIS_OBESE = list(
      description = "Obese status",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Ge 2020 Methods 'Covariate selection'); the paper does not state the obesity definition."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    n_observations = 87L,
    age_range = "Postnatal age 0.01-19.13 years (median 8.89 years); 20 infants (PNA <= 2 years), 9 children (2-12 years), 21 adolescents (> 12 years)",
    age_median = "8.89 years (postnatal age)",
    weight_range = "5th-95th percentile 2.6-98.1 kg",
    weight_median = "23.5 kg",
    sex_female_pct = 48,
    race_ethnicity = c(White = 70, `Black or African American` = 28, Unknown = 2),
    disease_state = "Pediatric patients receiving metoclopramide per standard of care; indications (Table S1): gastroesophageal reflux 18, gastroparesis 9, headache 8, nausea/vomiting without emetogenic chemotherapy 6, gastric motility 3, feeding intolerance 3, emesis 2, esophageal atresia 1.",
    dose_range = "Intravenous bolus median 0.1 (range 0.07-0.2) mg/kg in 20 patients (2 of them also on continuous infusion at 0.001 and 0.007 mg/kg/min); oral median 0.1 (0.04-0.15) mg/kg in 13 patients; 2 patients received both IV and oral doses and 15 received IV, oral and nasogastric/nasojejunal/gastrostomy doses.",
    regions = "United States (12 Pediatric Trials Network sites)",
    renal_function = "Schwartz creatinine clearance median 109.7 (5th-95th percentile 5.0-225.1) mL/min/1.73 m^2, available in 31 of 50 subjects",
    notes = "Pediatric Trials Network 'Pharmacokinetics of Understudied Drugs Administered to Children Per Standard of Care' trial (POPS; NCT01431326). Opportunistic sampling: 87 quantifiable concentrations (1 BLQ sample excluded; LLOQ 1 ng/mL), median 1.5 (range 1-6) samples per patient; 17 infants <= 120 days PNA contributed 34 samples. 3 of 50 patients (6%) had surgery within 24 h before a sampling dose. Demographics in Table 1 are median (5th-95th percentile) at the first recorded dose. NONMEM 7.4.1, FOCE-I. Eta shrinkage on CL 12.4%, epsilon shrinkage 19.2%."
  )

  ini({
    # Structural parameters (Ge 2020 Table 2 and Eqs. 5-10). Clearances and
    # volumes are per 70 kg.
    lka <- log(0.4); label("First-order absorption rate constant, all enteral routes (1/h)") # Ge 2020 Table 2: Ka = 0.4 1/h (RSE 20.4%; bootstrap 0.2-0.8); Eq. 5
    lcl <- log(19.6); label("Clearance at 70 kg (L/h)") # Ge 2020 Table 2: CL = 19.6 L/hour/70 kg (RSE 9.6%; bootstrap 12.5-21.8); Eq. 6
    lvc <- log(42.9); label("Central volume of distribution at 70 kg (L)") # Ge 2020 Table 2: Vc = 42.9 L/70 kg (RSE 25.2%; bootstrap 4.0-111.8); Eq. 7
    lq <- log(57.1); label("Intercompartmental clearance at 70 kg (L/h)") # Ge 2020 Table 2: Q = 57.1 L/hour/70 kg (RSE 24.0%; bootstrap 11.5-238.0); Eq. 8
    lvp <- log(83.9); label("Peripheral volume of distribution at 70 kg (L)") # Ge 2020 Table 2: Vp = 83.9 L/70 kg (RSE 12.6%; bootstrap 41.6-275.8); Eq. 9
    lfdepot <- log(0.97); label("Bioavailability, all enteral routes (fraction)") # Ge 2020 Table 2: F = 0.97 (RSE 12.5%; bootstrap 0.6-1.0); Eq. 10 'F = 97%'

    # Allometric exponents on WT/70, fixed a priori (Ge 2020 Methods
    # 'Covariate selection'; Results: 'Allometric exponents of 0.75 and 1 were
    # fixed for clearance (CL and Q) and volume of distribution (Vc and Vp)').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT on CL and Q (unitless)") # Ge 2020 Eqs. 6 and 8: (WT/70)^0.75, fixed
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT on Vc and Vp (unitless)") # Ge 2020 Eqs. 7 and 9: (WT/70)^1, fixed

    # IIV on CL only (IIV on Ka, Vc, Vp, F and Q gave eta shrinkage > 90% and
    # was dropped; Results). Reported as %CV; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.16501 # 42.4 %CV -> log(1 + 0.424^2) = 0.16501; Ge 2020 Table 2 'IIV, CL, %CV' = 42.4 (RSE 13.8%; bootstrap 25.7-57.4); eta shrinkage 12.4%

    # Residual error: proportional only (Results).
    propSd <- 0.333; label("Proportional residual error (fraction)") # Ge 2020 Table 2: proportional error = 33.3% (RSE 21.8%; bootstrap 23.6-39.7); epsilon shrinkage 19.2%
  })
  model({
    # Individual PK parameters (Ge 2020 Eqs. 1 and 5-10).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # A single bioavailability for oral, nasogastric, nasojejunal and
    # gastrostomy doses; intravenous doses go directly into central.
    f(depot) <- exp(lfdepot)

    # Dose in mg and volume in L give mg/L; x 1000 -> ng/mL, the units of the
    # assay (LLOQ 1 ng/mL) and of the paper's exposure targets.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
