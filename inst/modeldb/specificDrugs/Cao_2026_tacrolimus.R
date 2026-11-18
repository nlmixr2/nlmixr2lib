Cao_2026_tacrolimus <- function() {
  description <- "One-compartment population PK model for oral tacrolimus with first-order absorption in Chinese adults with membranous nephropathy, built from routine trough therapeutic-drug-monitoring concentrations (Cao 2026). Apparent oral clearance and apparent volume of distribution are allometrically scaled on body weight (fixed exponents 0.75 and 1, 70 kg reference); the absorption rate constant is fixed at 4.48 1/h from the literature because only troughs were sampled. None of the screened laboratory values or 18 concomitant medications was retained. Between-subject variability on CL/F only; additive residual error."
  reference <- paste(
    "Cao P, Jin YW, Duan YY, Wang J, Ping ZH, Zhao SY, Du JY, He SM, Zhang HB, Wang DD.",
    "Model-Informed Precision Dosing of Tacrolimus Initial Dosage Optimization in",
    "Membranous Nephropathy Patients. Drug Des Devel Ther. 2026;20:623781.",
    "doi:10.2147/DDDT.S623781.",
    "Final model Equations (6) and (7); parameter estimates Table 3.",
    "The fixed absorption rate constant is quoted from the paper's references 33 and 34:",
    "Wang DD, Chen X, Fu M, Zheng QS, Xu H, Li ZP. Model extrapolation to a real-world",
    "dataset: evaluation of tacrolimus population pharmacokinetics and drug interaction in",
    "pediatric liver transplantation patients. Xenobiotica. 2020;50(4):371-379.",
    "doi:10.1080/00498254.2019.1631505. And Wang D, Chen X, Xu H, Li Z. Population",
    "pharmacokinetics and dosing regimen optimization of tacrolimus in Chinese pediatric",
    "hematopoietic stem cell transplantation patients. Xenobiotica. 2020;50(2):178-185.",
    "doi:10.1080/00498254.2019.1601791.",
    sep = " "
  )
  vignette <- "Cao_2026_tacrolimus"
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
        "Equation 3, and Results, Equations 6-7). Cohort mean 69.65 kg (SD 11.50), median",
        "67.50, range 52.00-103.00 (Table 1). Treated as a baseline, time-fixed covariate.",
        "The paper extrapolates the adult model to children of 10-40 kg for its dosing",
        "simulations (Figures 6-7); no child contributed data to the fit."
      ),
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 49.35 years (SD 11.07), median 55.00, range 30.00-67.00 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained; 22 men / 9 women (Table 1, 'Gender (men/women)')."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 34.05 g/L (SD 6.79), median 34.50, range 19.00-44.50 (Table 1)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 23.04 IU/L (SD 14.46), median 17.00, range 7.00-65.00 (Table 1)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 21.38 IU/L (SD 7.71), median 20.00, range 13.00-43.00 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 71.64 umol/L (SD 25.28), median 68.00, range 38.00-154.00 (Table 1)."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 7.19 mmol/L (SD 3.02), median 6.61, range 3.17-14.98 (Table 1). The paper's Table 1 row is labelled 'Urea (mmol/L)'."
    ),
    TPRO = list(
      description = "Total protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 54.58 g/L (SD 8.34), median 56.10, range 38.00-66.50 (Table 1)."
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 5.87 mmol/L (SD 1.76), median 5.19, range 3.89-9.49 (Table 1)."
    ),
    TRIG = list(
      description = "Serum triglyceride",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 2.85 mmol/L (SD 2.57), median 2.38, range 1.05-13.87 (Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 3.15 umol/L (SD 0.92), median 3.00, range 1.60-4.50 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 8.21 umol/L (SD 3.22), median 7.10, range 3.90-15.60 (Table 1, row 'Total bilibrubin' (sic))."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 39.20% (SD 4.61), median 39.80, range 28.30-47.00 (Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained; cohort mean 128.90 g/L (SD 16.68), median 128.00, range 86.00-155.00 (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31,
    n_studies = 1,
    age_mean = "49.35 years (SD 11.07)",
    age_median = "55.00 years",
    age_range = "30.00-67.00 years",
    weight_mean = "69.65 kg (SD 11.50)",
    weight_median = "67.50 kg",
    weight_range = "52.00-103.00 kg",
    sex_female_pct = 29.0,
    disease_state = "Membranous nephropathy treated with oral tacrolimus.",
    dose_range = "Not reported. Concentrations came from routine therapeutic drug monitoring and the paper tabulates no administered tacrolimus doses or dosing frequencies.",
    regions = "China (single centre: Affiliated Hospital of Xuzhou Medical University, Jiangsu)",
    notes = paste(
      "Retrospective analysis of tacrolimus therapeutic-drug-monitoring concentrations",
      "collected between May 2022 and September 2025 (Methods, Enrolled patients).",
      "All observations are sparse TROUGH concentrations, sampled about 30 minutes before",
      "the next dose. Tacrolimus was assayed by enzyme-multiplied immunoassay (Emit 2000,",
      "Siemens; routine range 2-30 ng/mL).",
      "Baseline demographics and laboratory values are Table 1; concomitant medication",
      "counts are Table 2. The 18 concomitant medications screened were allisartan",
      "isoproxil, alprostadil, atorvastatin, Bailing capsule, co-trimoxazole, dapagliflozin,",
      "dipyridamole, entecavir, Fufang Shenyan tablet, Huangkui capsule, Jinshuibao tablet,",
      "levothyroxine, papaverine, prednisolone acetate, prednisone acetate, rivaroxaban,",
      "Shenkang injection and tripterygium glycosides; none was retained (Discussion).",
      "Globulin (20.80 g/L, SD 2.89), mean corpuscular hemoglobin (29.81 pg, SD 2.64) and",
      "mean corpuscular hemoglobin concentration (328.31 g/L, SD 10.94) were also screened",
      "and have no canonical covariate-column entry, so they are recorded only here.",
      "The model was qualified by goodness-of-fit plots, individual plots, a VPC",
      "(Figures 1-3) and a bootstrap (Table 3). The dosing simulations (Figures 4-7,",
      "Table 4) cover adults of 50-100 kg and, by extrapolation, children of 10-40 kg."
    )
  )

  ini({
    # Structural parameters -- Results, Equations (6) and (7); point estimates
    # also in Table 3. Reference subject: 70 kg.
    lcl <- log(9.54); label("Apparent oral clearance CL/F at 70 kg (L/h)")  # Table 3 / Eq. 6: CL/F = 9.54 L/h (bootstrap median 9.66 [1.13, 13.75])
    lvc <- log(76.3); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 3 / Eq. 7: V/F = 76.3 L (bootstrap median 78.9 [3.0, 214.0])

    # Ka was not estimated. Methods, Modeling: Ka was "fixed at 4.48 h-1 based
    # on references" 33 and 34 (Wang 2020, two Xenobiotica tacrolimus popPK
    # analyses); Table 3 prints "4.48 (fixed)" with no bootstrap row.
    lka <- fixed(log(4.48)); label("First-order oral absorption rate constant ka (1/h)")  # Table 3: Ka = 4.48 1/h (fixed)

    # Allometric exponents. Methods, Equation (3): "M was the allometric
    # coefficient: 0.75 for the CL/F and 1 for the V/F", citing Anderson &
    # Holford 2008. Asserted, not estimated -- neither appears in Table 3.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for WT/70 (unitless)")  # Methods Eq. 3 and Results Eq. 6 exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on V/F for WT/70 (unitless)")  # Methods Eq. 3 and Results Eq. 7 exponent 1

    # Between-subject variability. Methods Equation (1): B_i = TV(B) * exp(eta_i).
    # Table 3 labels the row 'omega_CL/F' (not omega^2); 0.172 is read as the
    # STANDARD DEVIATION. With that reading the model reproduces the paper's
    # target-attainment curves (Figures 5 and 7) to within a few percentage
    # points; reading 0.172 as a variance (SD 0.415) flattens every curve far
    # below the published 70-85% peaks. See the vignette.
    etalcl ~ 0.172^2  # Table 3: omega_CL/F = 0.172 (bootstrap median 0.094 [0.003, 0.381]), read as an SD

    # Residual error. Methods Equation (2): E_i = G_i + eps_1, purely additive.
    # Table 3 footnote: 'sigma_1, residual variability, additive error'. Read as
    # an SD in ng/mL, on the same reading as omega above.
    addSd <- 2.504; label("Additive residual error (ng/mL)")  # Table 3: sigma_1 = 2.504 (bootstrap median 2.601 [0.010, 3.508])
  })

  model({
    # Results Equations (6) and (7):
    #   CL/F = 9.54 * (weight/70)^0.75
    #   V/F  = 76.3 * (weight/70)
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

    Cc ~ add(addSd)
  })
}
