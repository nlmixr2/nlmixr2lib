Chen_2022b_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral tacrolimus whole-blood concentrations ",
    "in Chinese children with severe combined immunodeficiency (SCID) ",
    "undergoing haematopoietic stem cell transplantation (Chen 2022). The ",
    "absorption rate constant ka is fixed at 4.48 1/h from earlier paediatric ",
    "tacrolimus models. Apparent oral clearance CL/F is allometrically scaled ",
    "by body weight (fixed exponent 0.75, reference 70 kg) and apparent volume ",
    "V/F scales linearly with body weight (fixed exponent 1). Exponential IIV ",
    "on CL/F and V/F; combined proportional-plus-additive residual error."
  )
  reference <- paste0(
    "Chen X, Wang D, Zheng F, Zhai X, Xu H, Li Z. Population pharmacokinetics ",
    "and initial dose optimization of tacrolimus in children with severe ",
    "combined immunodeficiency undergoing hematopoietic stem cell ",
    "transplantation. Front Pharmacol. 2022;13:869939. ",
    "doi:10.3389/fphar.2022.869939"
  )
  vignette <- "Chen_2022b_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1), ",
        "both fixed, normalised to a standard weight of 70 kg (Chen 2022 ",
        "Methods Eq. 3; final-model Eqs. 6-7). Table 1: 7.28 +/- 1.62 kg, ",
        "median 7.50 (range 4.20-12.60) kg, so the 70 kg reference lies far ",
        "outside the observed range."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected (Table 1: 0.82 +/- 0.56 years, median 0.70, range 0.33-3.01) but not retained in the final model (Eqs. 6-7)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Collected (Table 1: 14 boys, 4 girls) but not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Collected (Table 1: median 33.20, range 25.10-40.80 g/L) but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Collected (Table 1: median 25.25, range 11.00-140.10 IU/L) but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Collected (Table 1: median 55.60, range 29.30-152.00 IU/L) but not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Collected (Table 1: median 18.00, range 9.00-27.00 umol/L) but not retained."
    ),
    BUN = list(
      description = "Urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Collected as urea (Table 1: median 2.45, range 0.60-7.00 mmol/L) but not retained."
    ),
    TPRO = list(
      description = "Total protein",
      units = "g/L",
      type = "continuous",
      notes = "Collected (Table 1: median 56.85, range 46.80-75.40 g/L) but not retained."
    ),
    TBA = list(
      description = "Total bile acid",
      units = "umol/L",
      type = "continuous",
      notes = "Collected (Table 1: median 5.90, range 0.10-21.30 umol/L) but not retained."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Collected (Table 1: median 2.40, range 0.80-24.40 umol/L) but not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Collected (Table 1: median 6.15, range 2.20-39.70 umol/L) but not retained."
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = "Collected (Table 1: median 26.61, range 22.80-53.31 %) but not retained, although haematocrit is a retained CL/F covariate in several other tacrolimus models."
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Collected (Table 1: median 87.10, range 69.00-167.00 g/L) but not retained."
    ),
    CONMED_STEROID = list(
      description = "Concomitant glucocorticoid",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1: 17 of 18 patients) but not retained (Discussion)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1: 13 of 18 patients) but not retained (Discussion)."
    ),
    CONMED_VANCOMYCIN = list(
      description = "Concomitant vancomycin",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1: 10 of 18 patients) but not retained (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    n_observations = 130L,
    age_range = "0.33-3.01 years",
    age_median = "0.70 years",
    weight_range = "4.20-12.60 kg",
    weight_median = "7.50 kg",
    sex_female_pct = 22.2,
    race_ethnicity = "Chinese (single centre, Children's Hospital of Fudan University, Shanghai).",
    disease_state = paste0(
      "Children with severe combined immunodeficiency (SCID) undergoing ",
      "haematopoietic stem cell transplantation and treated with oral ",
      "tacrolimus, February 2016 to April 2021, analysed retrospectively. ",
      "14 boys, 4 girls."
    ),
    dose_range = "Oral tacrolimus, doses adjusted by trough TDM; the individual doses are not tabulated.",
    regions = "China (Shanghai).",
    co_medication = paste0(
      "Caspofungin 9, ethambutol 10, glucocorticoids 17, isoniazid 14, ",
      "micafungin 9, mycophenolic acid 6, omeprazole 13, vancomycin 10 of 18 ",
      "(Table 1); none retained as a covariate."
    ),
    assay = "Whole-blood tacrolimus by Emit 2000 Tacrolimus Assay (Siemens), range 2.0-30 ng/mL.",
    notes = paste0(
      "Routine therapeutic-drug-monitoring data: 130 concentrations, mean ",
      "7.2 per patient. Table 1 also reports mean corpuscular haemoglobin ",
      "and its concentration, neither retained. NONMEM 7 FOCE-I; ",
      "1000-replicate bootstrap."
    )
  )

  ini({
    # Chen 2022 Table 2 'Parameter estimates of final model and bootstrap
    # validation'; final-model equations printed in Results as Eqs. 6-7:
    #   CL/F = 13.1 x (WT/70)^0.75
    #   V/F  = 10900 x (WT/70)
    lcl <- log(13.1); label("Apparent oral clearance CL/F at 70 kg (L/h)") # Table 2 CL/F = 13.1 L/h (SE 27.0%); Eq. 6
    # Table 2 prints V/F as 109 in units of 10^2 L; Eq. 7 prints 10900 L.
    lvc <- log(10900); label("Apparent volume of distribution V/F at 70 kg (L)") # Table 2 V/F = 109 x 10^2 L (SE 19.5%); Eq. 7
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)") # Table 2 Ka = 4.48 (fixed); Methods 'Ka was fixed at 4.48/h (Yang et al., 2015; Wang et al., 2019)'

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Methods Eq. 3, R = 0.75 for CL/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Methods Eq. 3, R = 1 for V/F

    # Table 2 omega values are read as SDs of eta. Replicating the Figure 4
    # target-attainment curves favours the SD reading over the variance
    # reading (see the vignette Assumptions section).
    etalcl ~ 0.203401 # Table 2 omega CL/F = 0.451 (SD), SE 44.8%
    etalvc ~ 0.350464 # Table 2 omega V/F = 0.592 (SD), SE 26.4%

    # Table 2 sigma 1 (proportional) and sigma 2 (additive) are read as SDs,
    # the same scale as omega in the same table; Methods Eq. 2
    # Mi = Ni x (1 + eps1) + eps2.
    propSd <- 0.257; label("Proportional residual error (fraction)") # Table 2 sigma 1 = 0.257 (SE 9.0%)
    addSd <- 1.265; label("Additive residual error (ng/mL)") # Table 2 sigma 2 = 1.265 (SE 17.5%)
  })

  model({
    ka <- exp(lka)
    # Eq. 6
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    # Eq. 7
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, vc in L -> mg/L; x 1000 gives ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
