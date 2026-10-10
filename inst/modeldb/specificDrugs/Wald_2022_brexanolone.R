Wald_2022_brexanolone <- function() {
  description <- paste(
    "Two-compartment intravenous population PK of brexanolone (exogenous",
    "allopregnanolone) in women with postpartum depression and healthy",
    "lactating women given the 60-h weight-based infusion, with allometric",
    "body-weight scaling (fixed exponents 0.75 and 1) and breast-milk",
    "concentration linked to plasma by an estimated milk-to-plasma",
    "concentration ratio."
  )
  reference <- paste(
    "Wald J, Henningsson A, Hanze E, Hoffmann E, Li H, Colquhoun H,",
    "Deligiannidis KM. Allopregnanolone Concentrations in Breast Milk and Plasma",
    "from Healthy Volunteers Receiving Brexanolone Injection, With Population",
    "Pharmacokinetic Modeling of Potential Relative Infant Dose. Clin",
    "Pharmacokinet. 2022;61(9):1307-1319. doi:10.1007/s40262-022-01155-w"
  )
  vignette <- "Wald_2022_brexanolone"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight; allometric scaling of CL, Q, V1 and V2",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference weight 82.9 kg is the median of the pooled N = 156 analysis",
        "population (Table 1; ESM Table S1 footnote). Fixed allometric",
        "exponents 0.75 for CL and Q and 1 for V1 and V2 (ESM 'PopPK Model",
        "Development'). The infusion rate is also weight-based (ug/kg/h), so",
        "the dose in ug is 'rate x WT x duration'."
      ),
      source_name = "WT"
    )
  )

  # Covariates screened by stepwise covariate modelling (ESM 'Final PopPK
  # Model') but not retained. Documentation only; not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "Screened on CL, V1 and V2; not retained (ESM 'Final PopPK Model')."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on CL, V1, V2 and Q; not retained (ESM 'Final PopPK Model')."
    ),
    RACE_BLACK = list(
      description = "Black or African American race (largest non-White group, 32.1%)",
      units = "(binary)",
      type = "binary",
      notes = "Race/ethnicity screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    ),
    CRCL = list(
      description = "Creatinine clearance (estimating equation not stated)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened on CL; not retained (ESM 'Final PopPK Model')."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "allopregnanolone",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "allopregnanolone",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 156L,
    n_studies = 5L,
    n_observations = paste(
      "1256 plasma + breast-milk concentrations in the final analysis from 155",
      "evaluable subjects (ESM 'PopPK Model Development')"
    ),
    age_range = "18-42 years",
    age_median = "27 years",
    weight_range = "44.9-150 kg",
    weight_median = "82.9 kg",
    sex_female_pct = 100,
    race_ethnicity = c(
      White = 50.6,
      Black = 32.1,
      Asian = 0.6,
      Other = 3.2,
      Hispanic = 13.5
    ),
    disease_state = paste(
      "Postpartum depression (144 patients from studies 202A, 202B, 202C and",
      "the open-label PPD study 201) and healthy lactating women <= 6 months",
      "postpartum (12 volunteers from study 547-CLP-108)"
    ),
    dose_range = paste(
      "Single 60-h continuous IV infusion titrated to 60 ug/kg/h (BRX60) or",
      "90 ug/kg/h (BRX90): 30 ug/kg/h 0-4 h, 60 ug/kg/h 4-24 h, 90 ug/kg/h",
      "24-52 h, 60 ug/kg/h 52-56 h, 30 ug/kg/h 56-60 h for BRX90"
    ),
    regions = "United States",
    notes = paste(
      "Baseline demographics in Table 1 (mean (SD) weight 85 (22) kg, BMI",
      "31.3 (7.8) kg/m^2). Plasma data from all five studies; milk data from",
      "547-CLP-108 only. Kp was estimated with the 547-CLP-108 plasma and milk",
      "data and held constant when the plasma model was updated with the full",
      "pooled dataset (ESM Table S1 footnote). Post hoc AUC and Cmax did not",
      "differ with vs without concomitant antidepressants (ESM Fig. S3).",
      "NONMEM 7.3."
    )
  )

  ini({
    lcl <- log(89.8)
    label("Clearance at 82.9 kg (L/h)")
    # ESM Table S1: CL = 89.8 L/h (RSE 1.8%; 95% CI 86.6-93.1)
    lvc <- log(117)
    label("Central volume of distribution at 82.9 kg (L)")
    # ESM Table S1: V1 = 117 L (RSE 22.9%; 95% CI 64.6-170)
    lq <- log(37.9)
    label("Intercompartmental clearance at 82.9 kg (L/h)")
    # ESM Table S1: Q = 37.9 L/h (RSE 7.7%; 95% CI 32.2-43.6)
    lvp <- log(470)
    label("Peripheral volume of distribution at 82.9 kg (L)")
    # ESM Table S1: V2 = 470 L (RSE 5.9%; 95% CI 415-524)
    lcmpr <- log(1.36)
    label("Milk-to-plasma concentration ratio (unitless)")
    # ESM Table S1: Kp (milk:plasma) = 1.36 (SE 0.235 on log scale; 95% CI 0.858-2.16)

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)")
    # ESM 'PopPK Model Development': allometric fixed constant 0.75 for CL and Q
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V1 and V2 (unitless)")
    # ESM 'PopPK Model Development': allometric fixed constant 1 for V1 and V2

    etalcl + etalvc ~ c(0.0435, 0.106, 1.15)
    # ESM Table S1: omega2 CL = 0.0435, omega2 CL,V1 = 0.106, omega2 V1 = 1.15 (variances; 21.1% and 147% CV)

    propSd <- 0.272
    label("Proportional residual error, plasma (fraction)")
    # ESM Table S1: proportional residual variability = 0.272
    propSd_Cmilk <- 0.272
    label("Proportional residual error, breast milk (fraction)")
    # ESM Table S1: single proportional residual variability = 0.272 for the joint plasma + milk fit
  })

  model({
    # ESM Table S1 footnote: P = theta * (WT / 82.9)^f, f = 0.75 for CL and Q,
    # f = 1 for V1 and V2.
    cl <- exp(lcl + etalcl) * (WT / 82.9)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 82.9)^e_wt_vc
    q <- exp(lq) * (WT / 82.9)^e_wt_cl
    vp <- exp(lvp) * (WT / 82.9)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in ug and volume in L give ug/L = ng/mL.
    Cc <- central / vc
    # Breast-milk concentration is the plasma concentration times the
    # estimated milk:plasma partition coefficient (main text 3.4; ESM Table S1).
    cmpr <- exp(lcmpr)
    Cmilk <- cmpr * Cc

    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
