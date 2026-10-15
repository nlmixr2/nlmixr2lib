Chen_2022_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral tacrolimus whole-blood concentrations in children with Crohn's ",
    "disease undergoing haematopoietic stem cell transplantation at a single ",
    "centre in China (Chen 2022). The absorption rate constant ka is fixed at ",
    "4.48 1/h from the literature. Apparent oral clearance CL/F is ",
    "allometrically scaled by body weight (fixed exponent 0.75, reference ",
    "70 kg) and reduced by 57% with concomitant posaconazole; apparent volume ",
    "V/F is scaled linearly by body weight (fixed exponent 1). Exponential ",
    "IIV on CL/F and V/F; combined proportional + additive residual error."
  )
  reference <- paste0(
    "Chen X, Wang D, Zheng F, Zhu L, Huang Y, Zhu Y, Huang Y, Xu H, Li Z. ",
    "Effects of posaconazole on tacrolimus population pharmacokinetics and ",
    "initial dose in children with Crohn's disease undergoing hematopoietic ",
    "stem cell transplantation. Front Pharmacol. 2022;13:758524. ",
    "doi:10.3389/fphar.2022.758524."
  )
  vignette <- "Chen_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to 70 kg (Chen 2022 Methods equation 3 and ",
        "Results equations 6-7). Table 1: 9.85 +/- 3.41 kg, median 9.50 ",
        "(range 3.70-20.60) kg."
      ),
      source_name = "WT"
    ),
    CONMED_POSACONAZOLE = list(
      description = "Concomitant posaconazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no posaconazole)",
      notes = paste0(
        "1 = the patient was treated with posaconazole, 0 = not (Chen 2022 ",
        "Results, text below equation 7). Fractional effect (1 + theta * POS) ",
        "on CL/F with theta = -0.57, i.e. CL/F ratio 1:0.43 without:with ",
        "posaconazole. Table 1: 12 of 51 patients received posaconazole. The ",
        "posaconazole formulation and dose are not reported, nor whether the ",
        "indicator was time-varying within a patient."
      ),
      source_name = "POS"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Methods, 'Covariate Model') but not retained in the final model. Table 1: 32 boys, 19 girls.",
      source_name = "Gender"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 1.86 +/- 1.38 years, median 1.36 (0.27-7.58).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 34.62 +/- 3.61 g/L, median 34.40 (27.50-43.30).",
      source_name = "Albumin"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1 (IU/L): 57.31 +/- 124.13, median 25.60 (7.70-789.40).",
      source_name = "Alanine transaminase"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1 (IU/L): 57.07 +/- 87.63, median 37.00 (12.50-628.20).",
      source_name = "Aspartate transaminase"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 17.55 +/- 3.80 umol/L, median 17.00 (11.00-28.00).",
      source_name = "Creatinine"
    ),
    UREA = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 2.92 +/- 1.27 mmol/L, median 2.90 (0.70-6.90). Reported as urea, not blood urea nitrogen.",
      source_name = "Urea"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 59.43 +/- 6.21 g/L, median 58.50 (48.10-72.20).",
      source_name = "Total protein"
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 7.47 +/- 7.18 umol/L, median 5.20 (0.90-33.90).",
      source_name = "Total bile acid"
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 3.45 +/- 6.99 umol/L, median 2.30 (0.80-51.80).",
      source_name = "Direct bilirubin"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 8.62 +/- 11.43 umol/L, median 6.60 (2.90-85.60).",
      source_name = "Total bilirubin"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 30.00 +/- 3.89%, median 30.10 (21.40-42.50).",
      source_name = "Hematocrit"
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 95.29 +/- 13.90 g/L, median 95.00 (66.00-147.00).",
      source_name = "Hemoglobin"
    ),
    MCH = list(
      description = "Mean corpuscular haemoglobin",
      units = "pg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 24.95 +/- 2.79 pg, median 25.10 (18.30-29.90).",
      source_name = "Mean corpuscular hemoglobin"
    ),
    MCHC = list(
      description = "Mean corpuscular haemoglobin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Table 1: 317.31 +/- 17.70 g/L, median 318.00 (280.00-348.00).",
      source_name = "Mean corpuscular hemoglobin concentration"
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic glucocorticoid indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no glucocorticoid)",
      notes = "Screened but not retained. Table 1: 40 of 51 patients.",
      source_name = "Glucocorticoids"
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mycophenolic acid)",
      notes = "Screened but not retained. Table 1: 26 of 51 patients.",
      source_name = "Mycophenolic acid"
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no omeprazole)",
      notes = "Screened but not retained. Table 1: 41 of 51 patients.",
      source_name = "Omeprazole"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 51L,
    n_studies = 1L,
    n_concentrations = 424L,
    age_range = "0.27-7.58 years (median 1.36; mean +/- SD 1.86 +/- 1.38)",
    weight_range = "3.70-20.60 kg (median 9.50; mean +/- SD 9.85 +/- 3.41)",
    sex_female_pct = 37.3,
    race_ethnicity = "Not reported (single centre in Shanghai, China)",
    disease_state = "Crohn's disease undergoing haematopoietic stem cell transplantation (HSCT), paediatric",
    dose_range = "Oral tacrolimus, initial dose 0.33-2 mg/day, then adjusted by TDM. Simulations used 0.1-0.8 mg/kg/day divided into two doses.",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    co_medication = "Posaconazole 12/51, glucocorticoids 40/51, omeprazole 41/51, mycophenolic acid 26/51 (Table 1).",
    sampling_design = paste0(
      "Retrospective routine TDM, October 2017 - December 2020; 424 whole-",
      "blood trough concentrations (about eight per patient) measured by ",
      "Emit 2000 immunoassay (range 2.0-30 ng/mL). Part of the clinical data ",
      "came from an earlier study by the same group (Wang 2020, Xenobiotica ",
      "50:178-85)."
    ),
    notes = paste0(
      "32 boys and 19 girls. Estimated in NONMEM 7 by FOCE-I; 1000-replicate ",
      "bootstrap and pcVPC. All concentrations were troughs, which is why a ",
      "one-compartment model with fixed ka was used."
    )
  )

  ini({
    # Chen 2022 Table 2 (final model), Results equations 6-7
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)")  # Methods 'PPK Model' and Table 2: Ka = 4.48 1/h (fixed), from Yang 2015 and Wang 2019
    lcl <- log(19.8); label("Apparent oral clearance CL/F at 70 kg without posaconazole (L/h)")  # Table 2 CL/F = 19.8 L/h, SE 7.0%; equation 6
    lvc <- log(11300); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 2 V/F = 113 (x10^2 L), SE 13.9%; equation 7 prints 11300

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods equation 3: power = 0.75 for CL/F (fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods equation 3: power = 1 for V/F (fixed)

    e_conmed_posaconazole_cl <- -0.57; label("Fractional change in CL/F with posaconazole (unitless)")  # Table 2 theta POS = -0.57, SE 12.1%; equation 6

    # Table 2 omega rows are the SD of eta (variance = omega^2). Their SEs
    # (15.8% and 14.7%) are below the sqrt(2/51) = 19.8% floor that any
    # variance estimate from 51 subjects must respect, and the SD reading
    # reproduces the paper's own Figure 5 and 6 simulations (see the
    # vignette).
    etalcl ~ 0.121801  # Table 2 omega CL/F = 0.349 (SD), SE 15.8%
    etalvc ~ 0.737881  # Table 2 omega V/F = 0.859 (SD), SE 14.7%

    propSd <- 0.259; label("Proportional residual error (fraction)")  # Table 2 sigma1 = 0.259, proportional error, SE 11.8%
    addSd <- 1.353; label("Additive residual error (ng/mL)")  # Table 2 sigma2 = 1.353, additive error, SE 13.2%
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_conmed_posaconazole_cl * CONMED_POSACONAZOLE)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    # Methods equation 2: O = IP * (1 + eps1) + eps2
    Cc ~ add(addSd) + prop(propSd)
  })
}
