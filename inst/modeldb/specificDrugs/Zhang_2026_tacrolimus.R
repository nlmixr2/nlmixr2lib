Zhang_2026_tacrolimus <- function() {
  description <- "One-compartment population PK model for oral tacrolimus with first-order absorption in Chinese adults with connective tissue disease-associated interstitial lung disease (CTD-ILD), built from routine trough therapeutic-drug-monitoring concentrations (Zhang 2026). Apparent oral clearance and apparent volume of distribution are allometrically scaled on body weight (fixed exponents 0.75 and 1, 70 kg reference); the absorption rate constant is fixed at 4.48 1/h from the literature because only troughs were sampled. No laboratory value or concomitant medication was retained (posaconazole was significant in screening but dropped as resting on one patient). Between-subject variability on CL/F only; combined proportional and additive residual error."
  reference <- paste(
    "Zhang GY, Jin YW, Wang J, Ding KW, Zhao SY, Ping ZH, Du JY, He SM, Zhang HB, Wang DD.",
    "Initial Dosage Optimization of Tacrolimus in Connective Tissue Disease-Associated",
    "Interstitial Lung Disease Patients: An Exploratory Study. Drug Des Devel Ther.",
    "2026;20:620915. doi:10.2147/DDDT.S620915.",
    "Final model Equations (8) and (9); parameter estimates Table 3.",
    "The fixed absorption rate constant is quoted from the paper's references 32-34:",
    "Wang DD, Chen X, Fu M, Zheng QS, Xu H, Li ZP. Xenobiotica. 2020;50(4):371-379.",
    "doi:10.1080/00498254.2019.1631505; Wang D, Chen X, Xu H, Li Z. Xenobiotica.",
    "2020;50(2):178-185. doi:10.1080/00498254.2019.1601791; Cai X, Song H, Jiao Z, et al.",
    "Eur J Pharm Sci. 2020;152:105448. doi:10.1016/j.ejps.2020.105448.",
    sep = " "
  )
  vignette <- "Zhang_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    # The paper names the assay (Emit 2000 Tacrolimus, Siemens), a whole-blood
    # immunoassay, but never states the matrix itself.
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling with a 70 kg standard weight, exponent 0.75 on CL/F and 1 on V/F,",
        "both fixed to the Anderson & Holford (2008) values rather than estimated (Methods,",
        "Equation 5, and Results, Equations 8-9). Cohort mean 66.07 kg (SD 8.94), median",
        "66.50, range 47.50-93.00 (Table 1). Treated as a baseline, time-fixed covariate."
      ),
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 3, dOFV -3.322); cohort mean 59.29 years (SD 10.23), median 55.00, range 42.00-75.00 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Table S1, model 2); 6 men / 16 women (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 4); cohort mean 40.33 g/L (SD 5.04), median 42.00, range 31.00-49.00 (Table 1)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 6); cohort mean 32.74 IU/L (SD 31.83), median 25.00, range 9.00-174.00 (Table 1)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 7); cohort mean 25.89 IU/L (SD 15.89), median 19.00, range 11.00-82.00 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 5); cohort mean 56.33 umol/L (SD 20.35), median 53.00, range 28.00-121.00 (Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 8); cohort mean 3.83 umol/L (SD 1.17), median 3.70, range 1.90-6.30 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 9); cohort mean 10.21 umol/L (SD 4.10), median 10.20, range 3.40-19.60 (Table 1, row 'Total bilibrubin' (sic))."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 10); cohort mean 44.07% (SD 4.71), median 43.60, range 34.30-53.70 (Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained (Table S1, model 11); cohort mean 144.16 g/L (SD 16.07), median 142.00, range 106.00-172.00 (Table 1)."
    ),
    CONMED_POSACONAZOLE = list(
      description = "Concomitant posaconazole indicator (1 = co-administered)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significant in forward inclusion (Table S1, model 25, dOFV -11.102) and retained by",
        "backward elimination (model 32, dOFV +11.102), but the authors left it out of the",
        "final model because the effect rested on one patient with two concentrations",
        "(Results, Modelling). The paper does not report the estimated effect size."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22,
    n_studies = 1,
    n_concentrations = 31,
    age_mean = "59.29 years (SD 10.23)",
    age_median = "55.00 years",
    age_range = "42.00-75.00 years",
    weight_mean = "66.07 kg (SD 8.94)",
    weight_median = "66.50 kg",
    weight_range = "47.50-93.00 kg",
    sex_female_pct = 72.7,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Connective tissue disease-associated interstitial lung disease treated with oral",
      "tacrolimus capsules: idiopathic inflammatory myopathy (n = 17), Sjogren syndrome",
      "(n = 3) and undifferentiated connective tissue disease (n = 2). Mean FVC 72.17%",
      "predicted (SD 14.40)."
    ),
    dose_range = "Not reported. Concentrations came from routine therapeutic drug monitoring and the paper tabulates no administered tacrolimus doses.",
    regions = "China (single centre: Affiliated Hospital of Xuzhou Medical University, Jiangsu)",
    notes = paste(
      "Retrospective analysis of tacrolimus therapeutic-drug-monitoring concentrations",
      "collected between August 2022 and September 2025 (Methods, Enrolled Patients).",
      "All 31 observations (1-3 per patient) are trough concentrations taken before the",
      "next dose, assayed by enzyme-multiplied immunoassay (Emit 2000, Siemens; routine",
      "range 2-30 ng/mL, LLOQ 2 ng/mL). One value (1.6 ng/mL) was below the LLOQ and was",
      "kept as reported; Table S2 shows deletion or LLOQ/2 imputation change CL/F by",
      "under 3%. Baseline demographics and laboratory values are Table 1; concomitant",
      "medication counts are Table 2 (15 drugs, most often prednisone 14/22 and",
      "co-trimoxazole 12/22). Table S1 records the 27 covariates screened on CL/F: the",
      "laboratory values in covariatesDataExcluded plus mean corpuscular hemoglobin",
      "(30.87 pg, SD 1.62) and its concentration (327.04 g/L, SD 8.15), which have no",
      "canonical covariate column, and all 15 concomitant medications. Globulin, urea,",
      "total protein, total cholesterol and triglyceride appear in Table 1 but not in the",
      "Table S1 screen. None of the 27 was retained. Prednisone",
      "reached p < 0.05 in the first forward step (Table S1, model 26) but not on top of",
      "posaconazole (model 30). The model was qualified by goodness-of-fit plots,",
      "individual plots, a VPC (Figures 1-3) and a bootstrap (Table 3, 99.6% successful).",
      "The dosing simulations (Figures 4-5, Table 4) cover adults of 50-100 kg."
    )
  )

  ini({
    # Structural parameters -- Results, Equations (8) and (9); point estimates
    # also in Table 3. Reference subject: 70 kg.
    lcl <- log(20.3); label("Apparent oral clearance CL/F at 70 kg (L/h)")  # Table 3 / Eq. 8: CL/F = 20.3 L/h, RSE 12.7% (bootstrap median 18.9 [12.2, 24.2])
    lvc <- log(1920); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 3 / Eq. 9: V/F = 1920 L, RSE 35.3% (bootstrap median 1905 [990, 6639])

    # Ka was not estimated: only troughs were sampled. Methods, Modeling: Ka
    # 'fixed at 4.48 h-1', citing references 32-34; Table 3 prints
    # '4.48 (fixed)' with no bootstrap row.
    lka <- fixed(log(4.48)); label("First-order oral absorption rate constant ka (1/h)")  # Table 3: Ka = 4.48 1/h (fixed)

    # Allometric exponents. Methods, Equation (5): 'F was the allometric
    # coefficient: 0.75 for the CL/F and 1 for the V/F', citing Anderson &
    # Holford 2008. Asserted, not estimated -- neither appears in Table 3.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for WT/70 (unitless)")  # Methods Eq. 5 and Results Eq. 8 exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on V/F for WT/70 (unitless)")  # Methods Eq. 5 and Results Eq. 9 exponent 1

    # Between-subject variability. Methods Equation (1): A_i = TV(A) * exp(eta_i).
    # Table 3 labels the row 'omega_CL/F' (not omega^2); 0.480 is read as the
    # STANDARD DEVIATION. With that reading the model reproduces the paper's
    # target-attainment curves (Figure 5) to a median 0.2 and at most 1.5
    # percentage points; reading 0.480 as a variance (SD 0.693) misses them by
    # a median 8 points. See the vignette.
    etalcl ~ 0.2304 # Table 3: omega_CL/F = 0.480 (RSE 20.8%) read as an SD; variance 0.480^2

    # Residual error. Methods Equation (4): B_i = C_i + C_i * eps_1 + eps_2.
    # Table 3 footnote: 'sigma1, residual variability, proportional error;
    # sigma2, residual variability, additive error'. Both read as SDs, on the
    # same reading as omega above.
    propSd <- 0.346; label("Proportional residual error (fraction)")  # Table 3: sigma1 = 0.346, RSE 31.8% (bootstrap median 0.339 [0.073, 0.516])
    addSd <- 0.322; label("Additive residual error (ng/mL)")  # Table 3: sigma2 = 0.322, RSE 100.5% (bootstrap median 0.411 [0.010, 1.456])
  })

  model({
    # Results Equations (8) and (9):
    #   CL/F = 20.3 * (weight/70)^0.75
    #   V/F  = 1920 * (weight/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment model with first-order oral absorption. Oral-only TDM
    # data, so cl and vc are the apparent quantities and no F is applied.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x1000 converts to ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
