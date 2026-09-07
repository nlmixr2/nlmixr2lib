Chen_2025_avibactam <- function() {
  description <- "One-compartment IV population PK model for the avibactam component of ceftazidime-avibactam in critically and non-critically ill Chinese adults with carbapenem-resistant Klebsiella pneumoniae infection (Chen 2025), with a median-normalized power-form creatinine-clearance effect on clearance."
  reference <- paste(
    "Chen Y, Chen B, Huang Y, Li X, Wu J, Lin R, Chen M, Liu M, Qiu H, Cheng Y.",
    "Population Pharmacokinetics-Based Evaluation of Ceftazidime-Avibactam",
    "Dosing Regimens in Critically and Non-Critically Ill Patients With",
    "Carbapenem-Resistant Klebsiella pneumoniae.",
    "Infect Drug Resist. 2025;18:941-953.",
    "doi:10.2147/IDR.S495279.",
    sep = " "
  )
  vignette <- "Chen_2025_ceftazidime_avibactam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated with the Cockcroft-Gault equation from serum creatinine recorded in the hospital electronic medical record (Chen 2025 Methods, Patients and Ethics)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw Cockcroft-Gault mL/min, NOT BSA-normalized. Stored under the canonical CRCL column per inst/references/covariate-columns.md, which accepts raw mL/min when the source paper applies no BSA normalization (same convention as Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R and the sibling combination-product pair Chandorkar_2015_ceftolozane.R / Chandorkar_2015_tazobactam.R). Reference value 71.3 mL/min = the cohort median (Chen 2025 Table 1; the paper states 71.3 is the median value of CrCL immediately below Eq. 3-4). Both analytes were fitted to the SAME 45 subjects, so unlike Chandorkar 2015 the two models share one centring constant. Cohort range 13.9-337.1 mL/min. Effect on CL is a median-normalized power form (CRCL / 71.3)^e_crcl_cl; the paper's typeset Eq. 3 places the fraction bar under CrCL^0.41 only, which cannot be the fitted form -- see the model file comment on e_crcl_cl and the validation vignette's Assumptions and deviations section. CrCL was the SOLE covariate retained after stepwise forward inclusion / backward elimination.",
      source_name        = "CrCL"
    )
  )

  # Screened during covariate selection but NOT retained in the final model
  # (Chen 2025 Methods, PopPK Modeling; Results, PopPK Modeling). Documented
  # here so the paper's covariate screen is preserved without declaring
  # covariates the model() block never references.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      reference_category = NULL,
      notes = "Median 62.0 kg (range 35.0-80.0), Chen 2025 Table 1. Screened; not retained. The final model carries no allometric term.",
      source_name = "Weight"
    ),
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      reference_category = NULL,
      notes = "Median 59.0 years (range 18.0-94.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator", units = "", type = "categorical",
      reference_category = "0 = male",
      notes = "9 of 45 subjects (20%) female, Chen 2025 Table 1 (reported as male/female 36/9). Screened; not retained.",
      source_name = "Gender"
    ),
    HT = list(
      description = "Body height", units = "cm", type = "continuous",
      reference_category = NULL,
      notes = "Median 170.0 cm (range 150.0-180.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "Height"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score", units = "", type = "continuous",
      reference_category = NULL,
      notes = "Used to define the disease-severity strata (> 15 critically ill, <= 15 non-critically ill; 25 vs 20 subjects, Chen 2025 Table 1). Screened; not retained -- Chen 2025 Discussion states explicitly that APACHE II scores do not significantly impact the PK of CAZ-AVI. The critically-ill / non-critically-ill split therefore drives only the PK/PD TARGET in the Monte Carlo simulation, not any structural or covariate term.",
      source_name = "APACHE II"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy status", units = "", type = "categorical",
      reference_category = "0 = no CRRT",
      notes = "15 of 45 subjects (33.3%) received CRRT, Chen 2025 Table 1. Screened; not retained. Chen 2025 Discussion attributes the null effect to the CRRT subgroup retaining substantial residual renal function (median CrCL 69.3 mL/min) and cites reports that CRRT alters clearance only below about 10 mL/min. Circuit settings when used: blood flow 160 mL/min, pre-replacement 2000 L/h, post-replacement 1000 L/h (Chen 2025 Results, Patients -- reproduced as printed).",
      source_name = "CRRT"
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 35.0 g/L (range 26.1-53.7), Chen 2025 Table 1. Screened; not retained.",
      source_name = "ALB"
    ),
    CREAT = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 79.0 umol/L (range 17.0-377.0), Chen 2025 Table 1. Screened as a covariate in its own right; the retained renal term is the Cockcroft-Gault CRCL derived from it.",
      source_name = "SCR"
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 13.4 umol/L (range 1.6-125.6), Chen 2025 Table 1. Screened; not retained.",
      source_name = "TBIL"
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 22.0 U/L (range 6.0-199.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 27.0 U/L (range 6.0-417.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "AST"
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 56.0 U/L (range 11.0-639.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "GGT"
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Median 113.0 U/L (range 26.0-523.0), Chen 2025 Table 1. Screened; not retained.",
      source_name = "ALP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = "18-94 years",
    age_median     = "59 years",
    weight_range   = "35.0-80.0 kg",
    weight_median  = "62.0 kg",
    sex_female_pct = 20,
    race_ethnicity = "Not reported by category; single-centre Chinese cohort (Fujian Medical University Union Hospital, Fuzhou)",
    disease_state  = "Adults with verified carbapenem-resistant Klebsiella pneumoniae infection receiving intravenous ceftazidime-avibactam. Infection sites: hospital-acquired pneumonia including ventilator-associated pneumonia 31 (68.9%), complicated intra-abdominal infection 5 (11.1%), bacteraemia associated with these infections 5 (11.1%), complicated urinary tract infection including pyelonephritis 4 (8.9%). Disease severity by APACHE II: 25 (55.6%) critically ill (> 15) and 20 (44.4%) non-critically ill (<= 15).",
    dose_range     = "Ceftazidime-avibactam 1.25 g or 2.5 g intravenously every 8, 12 or 24 h, each as a 2-hour infusion; the 4:1 fixed-ratio product delivers 250 mg or 500 mg of avibactam per dose. Most common regimen 2.5 g q8h (29 of 45 subjects).",
    regions        = "China (single centre, Fuzhou, Fujian)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 71.3 mL/min, range 13.9-337.1; 29 subjects (64.4%) with renal insufficiency (CrCL < 90 mL/min), 9 (20.0%) with CrCL 90-130, and 7 (15.6%) with augmented renal clearance (CrCL >= 130 mL/min). 15 subjects (33.3%) received continuous renal replacement therapy (median CrCL in that subgroup 69.3 mL/min).",
    notes          = "Prospective single-centre study, July 2021 to September 2023; baseline demographics per Chen 2025 Table 1. 91 steady-state plasma concentrations from 45 adults, of which 33 trough and 31 peak samples; sparse design of 1-3 samples per subject collected at or after the sixth dose. Assay LC-MS/MS validated to FDA bioanalytical guidance. Avibactam concentrations spanned 0.3-104.5 mg/L across sampling times 0-12 h post dose start (Chen 2025 Table 1, continued). Model estimated in Phoenix NLME 8.1 by first-order conditional estimation with extended least squares, separately from the ceftazidime model (distinct objective function values are reported for each analyte)."
  )

  ini({
    # Structural parameters for the reference subject (CRCL = 71.3 mL/min, the
    # cohort median). Chen 2025 Table 3 final-model column and Eq. 3-4.
    lcl <- log(3.09);  label("Clearance at CRCL 71.3 mL/min (L/h)")  # Chen 2025 Table 3: tvCL = 3.09 L/h (RSE 12.79%; bootstrap median 3.09, 95% CI 2.38-4.00); also Eq. 3
    lvc <- log(18.25); label("Volume of distribution (L)")           # Chen 2025 Table 3: tvV = 18.25 L (RSE 13.80%; bootstrap median 18.51, 95% CI 14.40-24.33); also Eq. 4

    # Median-normalized power effect of creatinine clearance on CL:
    #   CL = tvCL * (CRCL / 71.3)^e_crcl_cl * exp(etalcl)
    # NOT wrapped in fixed(): the table's "dCLdCrCL" row is labelled a
    # "fixed-parameter coefficient", but in Phoenix NLME "fixed parameter"
    # means a FIXED EFFECT (a THETA) as opposed to a random effect. The row
    # carries an RSE of 29.08% and a bootstrap 95% CI of 0.17-0.69, which
    # could only arise from estimation.
    #
    # Chen 2025 Eq. 3 is typeset as 3.09 x (CrCL^0.41 / 71.3) x exp(etaCL) --
    # the fraction bar sits under the exponentiated CrCL alone. That reading is
    # arithmetically impossible: at the median CrCL = 71.3 it gives
    # 3.09 x 71.3^0.41 / 71.3 = 0.25 L/h, an order of magnitude below the value
    # the same sentence calls "the typical value of CL". The median-normalized
    # power form used here returns exactly 3.09 L/h at CRCL = 71.3, matches the
    # Methods statement that "covariates were added to the model after median
    # normalization or as power functions", and reproduces the paper's own
    # external comparison (Discussion: a cohort averaging eGFR 50 mL/min
    # reported CL_AVI 4.9 L/h; this model gives 2.68 L/h at CRCL 50, the same
    # order of magnitude, whereas the literal reading gives 0.22 L/h).
    e_crcl_cl <- 0.41; label("CRCL exponent on CL (unitless)")  # Chen 2025 Table 3: dCLdCrCL = 0.41 (RSE 29.08%; bootstrap median 0.44, 95% CI 0.17-0.69); also Eq. 3

    # Inter-individual variability. Chen 2025 Table 3 headers the rows
    # "omega^2 CL" / "omega^2 V", but the tabulated numbers are STANDARD
    # DEVIATIONS on the log scale (Phoenix NLME reports the omega diagonal as
    # a stdev), which the paper itself converts to %CV by multiplying by 100:
    # Results, PopPK Modeling states that adding CrCL cut IIV in CL "from
    # 72.73% and 82.40% to 55.71% and 66.69%" for CAZ and AVI, and 0.5571 /
    # 0.6669 are exactly the Table 2 / Table 3 CL entries 0.56 / 0.67 before
    # rounding. No function of a VARIANCE of 0.67 yields 66.69% -- sqrt(0.67)
    # is 81.9% and sqrt(exp(0.67) - 1) is 97.7% -- so the header label is a
    # typesetting artefact. The same reading applies to the proportional
    # residual rows, whose 0.30 / 0.32 are 30% / 32% residual CV; read as
    # variances they would imply 55-57% residual error, far above the
    # LC-MS/MS assay imprecision the Methods describe.
    # nlmixr2 ini() takes VARIANCES, hence the squares below.
    etalcl ~ 0.4489 # 0.67^2; omega row for CL, Chen 2025 Table 3 (shrinkage 5.31%, RSE 17.38%) = 66.69% CV per Results text
    etalvc ~ 0.2601 # 0.51^2; omega row for V,  Chen 2025 Table 3 (shrinkage 16.39%, RSE 30.69%)

    # Residual variability. Chen 2025 Methods: "A proportional residual model
    # was applied to both models to assess residual variability."
    propSd <- 0.32; label("Proportional residual error (fraction)") # Chen 2025 Table 3: Proportional error = 0.32 (RSE 19.62%; bootstrap median 0.31, 95% CI 0.20-0.44)
  })
  model({
    # Individual PK parameters. CL carries the median-normalized power-form
    # creatinine-clearance effect (Chen 2025 Eq. 3); V has no covariate
    # (Chen 2025 Eq. 4).
    cl <- exp(lcl + etalcl) * (CRCL / 71.3)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
