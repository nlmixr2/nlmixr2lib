Chen_2020_ciclosporin <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral ciclosporin in Chinese children after bone marrow transplantation, with allometric body weight on CL/F and V/F and a median-normalised power effect of days post-transplant on CL/F (Chen 2020)"
  reference <- "Chen X, Yu X, Wang DD, Xu H, Li Z. Initial dosage optimization of ciclosporin in pediatric Chinese patients who underwent bone marrow transplants based on population pharmacokinetics. Exp Ther Med. 2020;20:401-408. doi:10.3892/etm.2020.8732"
  vignette <- "Chen_2020_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Ciclosporin whole-blood trough concentrations were measured with the Emit
  # 2000 Cyclosporin Specific assay (Chen 2020 Methods, 'Drug administration
  # and concentration detection'), so the modelled matrix is whole blood.
  compartmentData <- list(
    depot = list(analyte = "ciclosporin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ciclosporin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power scaling on both CL/F and V/F, normalised to a",
        "standard 70 kg adult (Chen 2020 Methods 'Covariate model',",
        "Equation C: POW = 0.75 for CL/F and 1 for V/F; Results Equations F",
        "and G). Cohort weight mean 8.40 +/- 3.28 kg, median 7.60 (range",
        "5.20-25.60) kg (Chen 2020 Table I), so every subject is",
        "extrapolated far below the 70 kg reference."
      ),
      source_name = "WT"
    ),
    POD = list(
      description = "Days post-transplant (bone marrow transplantation)",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying within subject. Enters CL/F as the median-normalised",
        "power (POD / 51.5)^0.749 (Chen 2020 Results Equation F; Table II",
        "theta-POD), with 51.5 days the cohort median (Table I: mean",
        "61.16 +/- 40.16, median 51.50, range 1.00-188.00 days). The",
        "positive exponent means apparent clearance rises with time after",
        "transplant. POD must be strictly positive: the power form sends",
        "CL/F to zero at POD = 0, and the observed data start at POD = 1."
      ),
      source_name = "POD"
    )
  )

  # Covariates that Chen 2020 screened in the stepwise covariate search but
  # did not retain (forward inclusion dOFV > 3.84, backward elimination
  # dOFV > 6.64; Chen 2020 Methods 'Statistical analysis'). Baseline
  # summaries are from Chen 2020 Table I.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "13 male / 5 female (Chen 2020 Table I); screened, not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Mean 1.60 +/- 1.15, median 1.22 (range 0.29-6.49) years (Chen 2020 Table I); screened, not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Mean 34.07 +/- 5.51, median 34.40 (range 1.20-45.20) g/L (Chen 2020 Table I); screened, not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Mean 51.26 +/- 65.97, median 30.00 (range 1.00-439.00) IU/L (Chen 2020 Table I); screened, not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Mean 56.00 +/- 57.58, median 36.30 (range 5.80-392.00) IU/L (Chen 2020 Table I); screened, not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 17.59 +/- 5.23, median 16.00 (range 8.00-38.00) umol/L (Chen 2020 Table I); screened, not retained."
    ),
    BUN = list(
      description = "Blood urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Mean 3.62 +/- 1.86, median 3.45 (range 0.60-11.10) mmol/L (Chen 2020 Table I, 'Urea'); screened, not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Mean 57.74 +/- 8.49, median 58.80 (range 38.40-77.80) g/L (Chen 2020 Table I); screened, not retained."
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 21.62 +/- 32.09, median 10.85 (range 0.40-201.20) umol/L (Chen 2020 Table I); screened, not retained."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 20.73 +/- 45.60, median 4.60 (range 0.10-301.90) umol/L (Chen 2020 Table I); screened, not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Mean 33.30 +/- 64.46, median 11.60 (range 1.30-384.20) umol/L (Chen 2020 Table I); screened, not retained."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Mean 30.64 +/- 6.76, median 30.40 (range 10.50-44.30) % (Chen 2020 Table I); screened, not retained."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Mean 101.16 +/- 21.86, median 101.00 (range 35.00-149.00) g/L (Chen 2020 Table I); screened, not retained."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic corticosteroid administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant systemic corticosteroid)",
      notes = "Glucocorticoids in 15 of 18 patients (Chen 2020 Table I); screened, not retained."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant omeprazole)",
      notes = "Omeprazole in 16 of 18 patients (Chen 2020 Table I); screened, not retained."
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes = "Phenobarbital in 2 of 18 patients (Chen 2020 Table I); screened, not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "0.29-6.49 years",
    age_median = "1.22 years (mean +/- SD 1.60 +/- 1.15)",
    weight_range = "5.20-25.60 kg",
    weight_median = "7.60 kg (mean +/- SD 8.40 +/- 3.28)",
    sex_female_pct = 27.8,
    race_ethnicity = "Chinese",
    disease_state = "Chinese children (<16 years) who underwent bone marrow transplantation (for severe aplastic anemia or leukemia) and received ciclosporin; recipients of other transplants (liver, kidney) were excluded",
    dose_range = "Initial 14-100 mg/day orally, subsequently adjusted on clinical efficacy, adverse events and trough TDM concentration",
    regions = "China (Children's Hospital of Fudan University, Shanghai)",
    notes = paste(
      "Retrospective therapeutic-drug-monitoring (TDM) dataset, September",
      "2016 to September 2019. Whole-blood trough concentrations only (drawn",
      "immediately before the next dose; Emit 2000 Cyclosporin Specific",
      "assay); the number of concentration records is not reported. Because",
      "only troughs were available the absorption rate constant was fixed to",
      "a literature value and CL/F and V/F are apparent (F-scaled)",
      "parameters. Covariates screened but not retained: sex, age, weight on",
      "other parameters, albumin, alanine transaminase, aspartate",
      "transaminase, creatinine, urea, total protein, total bile acid,",
      "direct bilirubin, total bilirubin, hematocrit, hemoglobin, mean",
      "corpuscular hemoglobin, mean corpuscular hemoglobin concentration, and",
      "concomitant glucocorticoids (15/18), mycophenolate mofetil (7/18),",
      "omeprazole (16/18), phenobarbital (2/18) and tacrolimus (2/18) (Chen",
      "2020 Methods 'Covariate model'; Table I). Mycophenolate mofetil,",
      "tacrolimus, mean corpuscular hemoglobin and mean corpuscular",
      "hemoglobin concentration have no canonical covariate column and are",
      "recorded here in prose only."
    )
  )

  ini({
    # Structural parameters; allometric reference body weight 70 kg
    lka <- fixed(log(0.68))
    label("Absorption rate constant (Ka, 1/h); from the literature, Chen 2020 references 16, 21 (Ni 2013) and 25 (Wang 2019) (Chen 2020 Methods 'Population pharmacokinetic modeling'; Table II)") # Table II 'Ka (h-1) 0.680', held constant
    lcl <- log(29.2)
    label("Apparent oral clearance at WT = 70 kg and POD = 51.5 days (CL/F, L/h); Chen 2020 Table II") # Table II 'CL/F (l/h) 29.200'
    lvc <- log(6550)
    label("Apparent volume of distribution at WT = 70 kg (V/F, L); Chen 2020 Table II") # Table II 'V/F (l) 6,550.000'

    # Allometric exponents fixed by convention (Chen 2020 Methods 'Covariate
    # model', Equation C, reference 26)
    e_wt_cl <- fixed(0.75)
    label("Allometric (WT/70) exponent on CL/F (unitless); Chen 2020 Methods Equation C") # Methods Equation C: POW = 0.75 for CL/F
    e_wt_vc <- fixed(1.0)
    label("Allometric (WT/70) exponent on V/F (unitless); Chen 2020 Methods Equation C") # Methods Equation C: POW = 1 for V/F

    # Covariate effect on CL/F
    e_pod_cl <- 0.749
    label("Power exponent of (POD/51.5) on CL/F (unitless); Chen 2020 Table II") # Table II 'theta POD 0.749'

    # Inter-individual variability. Chen 2020 Table II prints omega without
    # naming its scale. Read as standard deviations (the same convention the
    # same group's Wang 2019 cyclosporin paper in this journal confirms via
    # its quoted CV%s), because the variance reading implies a 97.5th
    # percentile near 1000 ng/mL at a typical 170 ng/mL, far outside the
    # published pcVPC (Chen 2020 Figure 2). See the vignette.
    etalcl ~ 0.393129
    # omega CL/F = 0.627, Chen 2020 Table II; 0.627^2 = 0.393129
    etalvc ~ 0.996004
    # omega V/F = 0.998, Chen 2020 Table II; 0.998^2 = 0.996004

    # Residual error: Y = F * (1 + eps1) + eps2 (Chen 2020 Methods 'Random
    # effect model', Equation B), combined proportional plus additive, on the
    # same standard-deviation scale as the omega column.
    propSd <- 0.447
    label("Proportional residual error (fraction); Chen 2020 Table II") # Table II 'sigma 1 0.447'
    addSd <- 70.071
    label("Additive residual error (ng/mL); Chen 2020 Table II") # Table II 'sigma 2 70.071'
  })

  model({
    # Individual PK parameters (Chen 2020 Results Equations F and G):
    #   CL/F = theta_CL/F * (WT/70)^0.75 * (POD/51.5)^theta_POD
    #   V/F  = theta_V/F  * (WT/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (POD / 51.5)^e_pod_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give central/vc in mg/L; x1000 converts to the
    # ng/mL whole-blood units the assay and Table III report.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
