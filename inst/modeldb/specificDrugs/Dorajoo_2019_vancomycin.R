Dorajoo_2019_vancomycin <- function() {
  description <- "One-compartment IV intermittent-infusion population PK model for vancomycin in adults with chronic kidney disease (CKD) not receiving renal replacement therapy (Dorajoo 2019). Clearance is a linear, cohort-mean-centred function of Cockcroft-Gault creatinine clearance (CL [L/h] = 1.30 x (1 + 0.023 x (CrCl [mL/min] - 33.8))) and volume of distribution is proportional to total body weight (V [L] = 1.23 x TBW [kg]). Exponential inter-individual variability on CL and V; combined additive plus proportional residual error. The final model was productized as the VancApp web-based dosing tool and prospectively implemented in routine care."
  reference <- "Dorajoo SR, Winata CL, Goh JHF, Ooi ST, Somani J, Yeoh LY, Lee SY, Yap CW, Chan A, Chae JW. Optimizing vancomycin dosing in chronic kidney disease by deriving and implementing a web-based tool using a population pharmacokinetics analysis. Front Pharmacol. 2019;10:641. doi:10.3389/fphar.2019.00641"
  vignette <- "Dorajoo_2019_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dorajoo 2019 Table 1 (construction cohort, n = 80): mean 57.8 kg, SD 15.7; Discussion reports median 55.8 kg (range 33.6-103.8). Source paper symbol TBW (total body weight); stored under the canonical WT column. The final model (Table 2) reports the volume of distribution in per-kilogram units (Vd = 1.23 L/kg), i.e. V [L] = 1.23 x TBW with no centering and no allometric exponent -- a purely multiplicative effect on V. The same linear-in-weight form appears in the authors' deposited VancApp server.R code (Supplementary R code: 'V <- 0.95*WT'), confirming the structure although that preliminary app build carries a different coefficient (see the vignette's Assumptions and deviations section).",
      source_name = "TBW"
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Dorajoo 2019 Methods ('Inclusion/Exclusion Criteria and Data Collection'): creatinine clearance estimated by the Cockcroft-Gault equation using TOTAL body weight ('to mirror actual clinical practice'). Stored under the canonical CRCL column with units mL/min (raw Cockcroft-Gault, not BSA-normalized), per the CRCL register entry's provision for raw mL/min when the source paper applies no BSA normalization, matching the Goti_2018_vancomycin.R, Buelga_2005_vancomycin.R and Moore_2016_vancomycin.R precedents. NOTE: Table 1 labels the row 'Creatinine clearance* (ml/min/1.73 m2)', but Cockcroft-Gault is NOT a BSA-normalized estimator; Table 3 labels the same quantity 'Creatinine clearance* (ml/min)' and the authors' deposited VancApp server.R code computes raw Cockcroft-Gault in mL/min with no BSA term. The '/1.73 m2' in Table 1 is therefore a labelling error and the column is raw mL/min. Construction-cohort mean 33.8 mL/min, SD 10.3 (Table 1). Cohort eligibility required CrCl < 60 mL/min, so the model's applicability domain is CrCl < 60. The centering constant 33.8 mL/min is the construction-cohort mean CrCl from Table 1; the paper does not print the centering value used in NONMEM (see the vignette's Assumptions and deviations section). Jelliffe-equation CrCl and MDRD / CKD-EPI eGFR were also screened as alternative renal-function descriptors (Supplementary Methods) but Cockcroft-Gault CrCl was retained in the final model.",
      source_name = "CrCl"
    )
  )

  # Covariates screened during stepwise covariate selection (Supplementary
  # Methods, 'List of covariates evaluated during model construction') that were
  # NOT retained in the final model. Results: 'CrCl significantly influenced the
  # clearance variability of vancomycin. All other potential covariates tested
  # were insignificant.' No point estimates are published for these, so they are
  # documentation only and are never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 71.7 years, SD 13.0 (Table 1); Discussion reports median 75 years (range 31-97)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Source paper reports 'gender'; screened on CL and Vd, not retained. Construction cohort 36.3% female (Table 1). Sex does enter the model's inputs indirectly through the Cockcroft-Gault CrCl calculation (x 0.85 for females), but not as a direct covariate on any PK parameter."
    ),
    HT = list(
      description = "Body height",
      units = "m",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 1.59 m, SD 0.08 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 22.6, SD 5.9 (Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 1.58, SD 0.23 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL and Vd as a standalone covariate; not retained (renal function entered the final model through Cockcroft-Gault CRCL instead). Construction cohort mean 150.6 umol/L, SD 73.6 (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 27.3 g/L, SD 6.7 (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 102.0 mg/L, SD 81.9 (Table 1)."
    ),
    PROCALCITONIN = list(
      description = "Serum procalcitonin",
      units = "ng/mL",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 10.1 ng/mL, SD 16.8 (Table 1)."
    ),
    WBC = list(
      description = "Total white blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened on CL and Vd; not retained. Construction cohort mean 13.76, SD 8.09 (Table 1). Neutrophil percentage (construction cohort mean 83.9%, SD 9.9; Table 1) was screened alongside it and likewise not retained; it has no canonical covariate column in inst/references/covariate-columns.md and is recorded here in prose rather than as a separate entry."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "median 75 years (range 31-97); mean 71.7, SD 13.0",
    age_median = "75 years",
    weight_range = "median 55.8 kg (range 33.6-103.8); mean 57.8, SD 15.7",
    weight_median = "55.8 kg",
    sex_female_pct = 36.3,
    race_ethnicity = "Not reported by category; the Discussion describes the study population as a 'multi-ethnic Asian cohort' at a single Singapore centre.",
    disease_state = "Adult inpatients with chronic kidney disease receiving intravenous vancomycin, with a baseline Cockcroft-Gault creatinine clearance < 60 mL/min and NOT receiving renal replacement therapy. 75% were at CKD stage 3 or beyond by MDRD eGFR (Table 1: CKD 1 6.3%, CKD 2 18.8%, CKD 3a 22.5%, CKD 3b 25.0%, CKD 4 23.8%, CKD 5 3.8%). Mean serum albumin 27.3 g/L and mean C-reactive protein 102.0 mg/L indicate a substantially inflamed, hypoalbuminaemic population.",
    dose_range = "Intravenous vancomycin by intermittent infusion under an in-house protocol: weight-based doses of 15-20 mg/kg at intervals of 12 or 24 h, infused at a maximum rate of 500 mg/h, with a maximum of 2 g per infusion regardless of total body weight.",
    regions = "Singapore (Khoo Teck Puat Hospital, a 600-bed tertiary institution)",
    renal_function = "Cockcroft-Gault CrCl (total body weight) mean 33.8 mL/min, SD 10.3; MDRD eGFR mean 46.2 mL/min/1.73 m^2, SD 23.0; serum creatinine mean 150.6 umol/L, SD 73.6 (Table 1). Patients on renal replacement therapy were excluded.",
    n_concentrations = 170L,
    notes = "Demographics from Dorajoo 2019 Table 1 (construction cohort: n = 80 patients, 170 vancomycin concentrations collected over the first 120 h of therapy, mean 2.1 +/- 1.3 concentrations per patient, range 1-4). Retrospective records of patients who received IV vancomycin from 1 April 2013 to 31 March 2014. Inclusion required at least two IV vancomycin doses over a 72-h period, a baseline Cockcroft-Gault CrCl < 60 mL/min, and at least one measured vancomycin concentration; exclusions were renal replacement therapy, missing baseline weight or serum creatinine, and IV vancomycin within the preceding 2 weeks. Modelling used NONMEM 7.3 with FOCE-I; 1000 bootstrap replicates supported the final estimates (Table 2 bootstrap column). A one-compartment model was preferred over two-compartment on AIC (744.5 vs 764.6) and on lower RSEs despite similar OFVs (704.5 vs 700.6; Supplementary Methods). The authors note that sampling was predominantly trough-based, which could have led to misspecification of Vd and a preferential fit of the one-compartment model. A SEPARATE temporal validation cohort (n = 112 patients, 289 concentrations, 1 April 2014 - 31 March 2015; Table 1) was used for external validation only and did NOT contribute to parameter estimation: mean absolute error 1.25 mg/L, mean absolute prediction error 8.6% (Supplementary Table S1). A further 33-patient post-implementation cohort (78 troughs) gave MAE 3.65 mg/L and MAPE 32.9% (Table 4)."
  )

  ini({
    # Structural parameters (Dorajoo 2019 Table 2, 'Final estimate' column).
    # The final model is:
    #   CL (L/h) = 1.30 x (1 + 0.023 x (CrCl [mL/min] - 33.8))
    #   V  (L)   = 1.23 x TBW [kg]
    # The linear, centred, fractional form of the CrCl effect is the form used in
    # the authors' own deposited VancApp server.R code (Supplementary R code:
    # 'CL <- 1.4*( 1 + 0.0224*(CRCL - 35.75))'); see the vignette's Assumptions
    # and deviations section for why the Table 2 values, not that preliminary
    # app build's values, are the ones packaged here, and for the digitisation
    # of Supplementary Figures S1 and S2 that confirms the form and the 33.8
    # mL/min centering constant.
    lcl <- log(1.30)
    label("Clearance at the cohort-mean creatinine clearance of 33.8 mL/min (L/h)") # Dorajoo 2019 Table 2: Cl = 1.30 L/h (RSE 7.2%, bootstrap 1.30, 95% CI 1.18-1.60)
    lvc <- log(1.23)
    label("Volume of distribution per kg total body weight (L/kg)") # Dorajoo 2019 Table 2: Vd = 1.23 L/kg (RSE 4.9%, bootstrap 1.23, 95% CI 1.12-1.37)

    # Covariate effect. Fractional change in CL per mL/min of Cockcroft-Gault
    # CrCl away from the 33.8 mL/min cohort mean. Over the cohort's eligible
    # CrCl range (0 to < 60 mL/min) the multiplier 1 + 0.023 x (CrCl - 33.8)
    # stays strictly positive (0.223 at CrCl = 0; 1.603 at CrCl = 60), so no
    # lower guard is required.
    e_crcl_cl <- 0.023
    label("Linear coefficient of creatinine clearance on CL (fractional change per mL/min)") # Dorajoo 2019 Table 2: theta CrCl = 0.023 (RSE 24.2%, bootstrap 0.023, 95% CI 0.011-0.033)

    # Inter-individual variability. Methods: 'The interindividual variability
    # (IIV) in the model was determined by exponential random effects.' Table 2
    # reports the etas as percentages, read here as coefficients of variation;
    # for a log-normal eta, omega^2 = log(CV^2 + 1):
    #   etalcl: 54.4% CV -> omega^2 = log(0.544^2 + 1) = log(1.295936) = 0.25925
    #   etalvc: 22.5% CV -> omega^2 = log(0.225^2 + 1) = log(1.050625) = 0.04939
    # Single-quoted source-table text below: rxode2 promotes an unlabelled
    # trailing ini() comment into label(), and a double quote would break it.
    etalcl ~ 0.25925 # Dorajoo 2019 Table 2, row 'eta Cl (%)' = 54.4 (RSE 10.2%, shrinkage 12%, bootstrap 53.8, 95% CI 41.2-65.5)
    etalvc ~ 0.04939 # Dorajoo 2019 Table 2, row 'eta Vd (%)' = 22.5 (RSE 17.1%, shrinkage 26%, bootstrap 21.9, 95% CI 12.6-29.3)

    # Residual error. Methods: 'the residual variability was modelled as a
    # mixture of additive and proportional error structures.' Table 2 reports
    # the additive term in mg/L (i.e. already a standard deviation) and the
    # proportional term as a variance (sigma^2), so propSd = sqrt(0.001).
    addSd <- 2.46
    label("Additive residual error (mg/L)") # Dorajoo 2019 Table 2: additive error = 2.46 mg/L (RSE 12.3%, bootstrap 2.44, 95% CI 1.85-2.98)
    propSd <- 0.0316228
    label("Proportional residual error (fraction)") # Dorajoo 2019 Table 2: proportional error sigma^2 = 0.001, so sigma = sqrt(0.001) = 0.0316
  })
  model({
    # Individual PK parameters. CRCL is the raw Cockcroft-Gault creatinine
    # clearance in mL/min, centred at the 33.8 mL/min construction-cohort mean;
    # WT is total body weight in kg and enters V linearly (Vd published in L/kg).
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 33.8))
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
