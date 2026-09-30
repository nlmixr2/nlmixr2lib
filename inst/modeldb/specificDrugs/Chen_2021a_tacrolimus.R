Chen_2021a_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral tacrolimus whole-blood concentrations in children with chronic ",
    "granulomatous disease (CGD) undergoing haematopoietic stem cell ",
    "transplantation at a single centre in China (Chen 2021). The absorption rate constant ka is fixed at ",
    "4.48 1/h from the literature. Apparent oral clearance CL/F is ",
    "allometrically scaled by body weight (fixed exponent 0.75, reference ",
    "70 kg) and reduced by 61.2% with concomitant voriconazole; apparent ",
    "volume V/F is scaled linearly by body weight (fixed exponent 1). ",
    "Exponential IIV on CL/F and V/F; combined proportional + additive ",
    "residual error."
  )
  reference <- paste0(
    "Chen X, Wang D, Lan J, Wang G, Zhu L, Xu X, Zhai X, Xu H, Li Z. Effects ",
    "of voriconazole on population pharmacokinetics and optimization of the ",
    "initial dose of tacrolimus in children with chronic granulomatous ",
    "disease undergoing hematopoietic stem cell transplantation. Ann Transl ",
    "Med. 2021;9(18):1477. doi:10.21037/atm-21-4124."
  )
  vignette <- "Chen_2021a_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to 70 kg (Chen 2021 Methods equation 3 and ",
        "Results equations 6-7). Table 1: 11.17 +/- 3.77 kg, median 10.00 ",
        "(range 6.30-24.80) kg."
      ),
      source_name = "WT"
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole)",
      notes = paste0(
        "1 = the patient took voriconazole, 0 = not (Chen 2021 Results, text ",
        "below equation 7). Fractional effect (1 + theta * VRC) on CL/F with ",
        "theta = -0.612, i.e. CL/F ratio 1:0.388 without:with voriconazole. ",
        "Table 1: 32 of 34 patients received voriconazole. The paper does not ",
        "state whether the indicator was time-varying within a patient."
      ),
      source_name = "VRC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Collected (Methods, Data collection) but not retained in the final model. Table 1: 33 boys, 1 girl.",
      source_name = "Gender"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 2.29 +/- 1.89 years, median 1.41 (0.38-9.28).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 36.99 +/- 2.97 g/L, median 37.20 (27.70-43.10).",
      source_name = "Albumin"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1 (IU/L): 25.36 +/- 15.38, median 21.40 (5.20-70.30).",
      source_name = "Alanine transaminase"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1 (IU/L): 41.70 +/- 37.39, median 31.00 (18.90-228.50).",
      source_name = "Aspartate transaminase"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 22.09 +/- 9.57 umol/L, median 20.00 (14.00-69.00).",
      source_name = "Creatinine"
    ),
    UREA = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 3.27 +/- 1.43 mmol/L, median 3.15 (1.00-8.50). Reported as urea, not blood urea nitrogen.",
      source_name = "Urea"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 60.91 +/- 5.68 g/L, median 61.15 (45.50-74.20).",
      source_name = "Total protein"
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 7.34 +/- 6.54 umol/L, median 5.10 (1.60-31.70).",
      source_name = "Total bile acid"
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 2.60 +/- 0.97 umol/L, median 2.50 (1.00-5.60).",
      source_name = "Direct bilirubin"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 7.94 +/- 3.16 umol/L, median 7.20 (2.70-17.60).",
      source_name = "Total bilirubin"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 26.47 +/- 4.08%, median 25.67 (16.60-38.10).",
      source_name = "Hematocrit"
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 88.74 +/- 13.80 g/L, median 85.50 (57.00-126.00).",
      source_name = "Hemoglobin"
    ),
    MCH = list(
      description = "Mean corpuscular haemoglobin",
      units = "pg",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 25.61 +/- 2.16 pg, median 25.50 (21.20-31.30).",
      source_name = "Mean corpuscular hemoglobin"
    ),
    MCHC = list(
      description = "Mean corpuscular haemoglobin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Collected but not retained. Table 1: 335.24 +/- 12.67 g/L, median 333.00 (315.00-374.00).",
      source_name = "Mean corpuscular hemoglobin concentration"
    ),
    CONMED_CASPOFUNGIN = list(
      description = "Concomitant caspofungin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no caspofungin)",
      notes = "Collected but not retained. Table 1: 21 of 34 patients.",
      source_name = "Caspofungin"
    ),
    CONMED_ETHAMBUTOL = list(
      description = "Concomitant ethambutol indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ethambutol)",
      notes = "Collected but not retained. Table 1: 26 of 34 patients.",
      source_name = "Ethambutol"
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic glucocorticoid indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no glucocorticoid)",
      notes = "Collected but not retained. Table 1: 23 of 34 patients.",
      source_name = "Glucocorticoids"
    ),
    CONMED_ISONIAZID = list(
      description = "Concomitant isoniazid indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no isoniazid)",
      notes = "Collected but not retained. Table 1: 31 of 34 patients.",
      source_name = "Isoniazid"
    ),
    CONMED_MICAFUNGIN = list(
      description = "Concomitant micafungin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no micafungin)",
      notes = "Collected but not retained. Table 1: 7 of 34 patients.",
      source_name = "Micafungin"
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mycophenolic acid)",
      notes = "Collected but not retained. Table 1: 10 of 34 patients.",
      source_name = "Mycophenolic acid"
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no omeprazole)",
      notes = "Collected but not retained. Table 1: 34 of 34 patients (no contrast available).",
      source_name = "Omeprazole"
    ),
    CONMED_VANCOMYCIN = list(
      description = "Concomitant vancomycin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no vancomycin)",
      notes = "Collected but not retained. Table 1: 14 of 34 patients.",
      source_name = "Vancomycin"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_concentrations = 293L,
    age_range = "0.38-9.28 years (median 1.41; mean +/- SD 2.29 +/- 1.89)",
    weight_range = "6.30-24.80 kg (median 10.00; mean +/- SD 11.17 +/- 3.77)",
    sex_female_pct = 2.9,
    race_ethnicity = "Not reported (single centre in Shanghai, China)",
    disease_state = "Chronic granulomatous disease (CGD) undergoing haematopoietic stem cell transplantation (HSCT), paediatric",
    dose_range = "Oral tacrolimus; administered dose range not reported. Simulations used 0.1-0.8 mg/kg/day divided into two doses.",
    regions = "China (single centre: Children's Hospital of Fudan University, Shanghai)",
    sampling_design = paste0(
      "Retrospective routine TDM, May 2016 - January 2021; 293 whole-blood ",
      "concentrations (mean 8.6 per patient) measured by Emit 2000 ",
      "immunoassay (range 2.0-30 ng/mL). Part of the data came from an ",
      "earlier study by the same group (Wang 2020, Xenobiotica 50:178-85)."
    ),
    notes = paste0(
      "33 boys and 1 girl. Estimated in NONMEM 7 by FOCE-I; 1000-replicate ",
      "bootstrap and pcVPC. Co-medications (Table 1): voriconazole 32/34, ",
      "omeprazole 34/34, isoniazid 31/34, ethambutol 26/34, glucocorticoids ",
      "23/34, caspofungin 21/34, vancomycin 14/34, mycophenolic acid 10/34, ",
      "micafungin 7/34."
    )
  )

  ini({
    # Chen 2021 Table 2 (final model), Results equations 6-7
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)")  # Methods 'Population pharmacokinetic model' and Table 2: Ka = 4.48 1/h (fixed), from refs 24-25
    lcl <- log(35.4); label("Apparent oral clearance CL/F at 70 kg without voriconazole (L/h)")  # Table 2 CL/F = 35.4 L/h, SE 12.2%; equation 6
    lvc <- log(5970); label("Apparent volume of distribution V/F at 70 kg (L)")  # Table 2 V/F = 59.7 (x10^2 L), SE 25.1%; equation 7 prints 5970

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods equation 3: power = 0.75 for CL/F (fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)")  # Methods equation 3: power = 1 for V/F (fixed)

    e_conmed_voriconazole_cl <- -0.612; label("Fractional change in CL/F with voriconazole (unitless)")  # Table 2 theta VRC = -0.612, SE 12.9%; equation 6

    # Table 2 omega rows are the SD of eta (variance = omega^2). Their SEs
    # (15.2% and 14.9%) are below the sqrt(2/34) = 24.3% floor that any
    # variance estimate from 34 subjects must respect, and the SD reading
    # reproduces the paper's own Figure 5 target-attainment curves (see the
    # vignette).
    etalcl ~ 0.275625  # Table 2 omega CL/F = 0.525 (SD), SE 15.2%
    etalvc ~ 0.680625  # Table 2 omega V/F = 0.825 (SD), SE 14.9%

    propSd <- 0.386; label("Proportional residual error (fraction)")  # Table 2 sigma1 = 0.386, proportional error, SE 10.1%
    addSd <- 0.354; label("Additive residual error (ng/mL)")  # Table 2 sigma2 = 0.354, additive error, SE 116.4%
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_conmed_voriconazole_cl * CONMED_VORICONAZOLE)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg, V/F L -> mg/L; x1000 to ng/mL
    Cc <- 1000 * central / vc
    # Methods equation 2: O = IPC * (1 + eps1) + eps2
    Cc ~ add(addSd) + prop(propSd)
  })
}
