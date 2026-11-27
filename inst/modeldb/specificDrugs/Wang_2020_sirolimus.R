Wang_2020_sirolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral sirolimus whole-blood trough concentrations in Chinese children ",
    "with tuberous sclerosis complex (Wang 2020). The absorption rate ",
    "constant ka is fixed at 0.485 1/h from the literature. Apparent oral ",
    "clearance CL/F and apparent volume V/F are allometrically scaled by ",
    "body weight (fixed exponents 0.75 and 1, reference 70 kg); no other ",
    "covariate, including the CYP3A4 and CYP3A5 genotypes, was retained. ",
    "Exponential IIV on CL/F and V/F; combined proportional and additive ",
    "residual error."
  )
  reference <- paste0(
    "Wang DD, Chen X, Xu H, Li ZP. Initial Dosage Recommendation for ",
    "Sirolimus in Children With Tuberous Sclerosis Complex. Front Pharmacol. ",
    "2020;11:890. doi:10.3389/fphar.2020.00890."
  )
  vignette <- "Wang_2020_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to a standard weight of 70 kg (Wang 2020 ",
        "Methods Eq. 3; Results Eqs. 7-8). Table 1: 23.83 +/- 8.88 kg, ",
        "median 22.00 (range 10.00-50.00) kg."
      ),
      source_name = "weight"
    )
  )

  # Covariates collected (Table 1, Table 2) and explored (Abstract; Methods
  # Eqs. 4-6) but not retained in the final model (Results Eqs. 7-8).
  # Documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Not retained. Table 1: 6.16 +/- 2.80, median 5.75 (1.08-13.95) years."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Not retained. Table 1: 7 boys, 8 girls."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 46.10 (40.30-48.30) g/L."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 7.00 (1.00-28.00) IU/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 24.00 (10.00-37.00) IU/L."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 29.00 (16.00-59.00) umol/L."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Not retained. Table 1 'Urea': median 4.50 (1.80-8.20) mmol/L."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 70.60 (64.70-81.10) g/L."
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 2.60 (1.30-11.80) umol/L."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 1.70 (0.40-4.80) umol/L."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 5.00 (1.00-17.70) umol/L."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Not retained. Table 1: median 37.75 (30.80-45.30) %."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Not retained. Table 1: median 128.50 (107.00-147.00) g/L."
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (1 = at least one *1 allele)",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Not retained (Discussion: the CYP3A4 and CYP3A5 polymorphisms 'were ",
        "not successfully included as covariates'). Table 2: CYP3A5 *1/*3 ",
        "in 7 and *3/*3 in 8 of 15 patients; CYP3A4 *1A/*1G in 6 of 15. ",
        "Eleven other genes (ABCB1, ABCC4, ABCC8, ABCG2, CYP2C9, CYP2C19, ",
        "CYP4F2, UGT1A1, UGT1A8, UGT2B15) were also genotyped (Table 2)."
      )
    )
  )

  # Methods 'TDM and Pharmacogenomic Analysis': Emit 2000 Sirolimus Assay
  # (Siemens), a whole-blood immunoassay, linear range 3.5-30 ng/mL.
  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "1.08-13.95 years (median 5.75; mean +/- SD 6.16 +/- 2.80)",
    weight_range = "10.00-50.00 kg (median 22.00; mean +/- SD 23.83 +/- 8.88)",
    sex_female_pct = 53.3,
    race_ethnicity = c(Asian = 100),
    disease_state = "Tuberous sclerosis complex (TSC), pediatric",
    dose_range = "Oral sirolimus titrated by therapeutic drug monitoring; individual doses not tabulated",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    sampling_design = paste0(
      "Retrospective routine TDM, June 2016 - September 2019; whole-blood ",
      "concentrations by Emit 2000 immunoassay (linear range 3.5-30 ng/mL). ",
      "The number of concentrations is not reported."
    ),
    notes = paste0(
      "Patients on CYP3A4 or P-gp inhibitors or inducers were excluded. ",
      "Part of the clinical data came from the group's earlier TSC analysis ",
      "(Wang 2018, Int J Clin Exp Med 11:12302). Estimated in NONMEM 7 by ",
      "FOCE-I; bootstrap (1,000 replicates) and pcVPC."
    )
  )

  ini({
    # Wang 2020 Table 3 (final model), Results Eqs. 7-8
    lka <- fixed(log(0.485)); label("Absorption rate constant ka (1/h)")  # Table 3 Ka = 0.485 (fixed); Methods, Population Pharmacokinetic Model
    lcl <- log(6.48); label("Apparent oral clearance CL/F at 70 kg (L/h)")  # Table 3 CL/F = 6.48 L/h, SE 33.2%; Eq. 7
    lvc <- log(124); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 3 V/F = 124 L, SE 69.8%; Eq. 8

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods Eq. 3, index 0.75 for CL/F; Eq. 7
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods Eq. 3, index 1 for V/F; Eq. 8

    # Table 3 omega rows are read as SDs of eta and squared here. The SD
    # reading reproduces the probability-of-target curves of Figures 2-3; the
    # variance reading flattens them (see the vignette).
    etalcl ~ 0.004356 # Table 3 omega CL/F = 0.066 (SD), SE 70.9%; 0.066^2
    etalvc ~ 0.000009 # Table 3 omega V/F = 0.003 (SD), SE 60.4%; 0.003^2

    # Methods Eq. 2: OBS = PRE * (1 + eps1) + eps2. Sigma rows read as SDs,
    # by the same notation as the omega rows.
    propSd <- 0.312; label("Proportional residual error (fraction)")  # Table 3 sigma1 = 0.312, proportional error, SE 27.2%
    addSd <- 1.249; label("Additive residual error (ng/mL)")  # Table 3 sigma2 = 1.249, additive error, SE 54.8%
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
