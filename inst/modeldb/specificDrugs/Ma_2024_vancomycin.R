Ma_2024_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in 245 elderly (>=65 years) inpatients undergoing therapeutic drug monitoring at a single center in Chongqing, China (Ma 2024). Clearance scales by power exponent with Cockcroft-Gault creatinine clearance (raw mL/min, reference 65.24); volume of distribution has no retained covariate. This is the popPK layer of a study whose headline product is a machine-learning ensemble that consumes the empirical-Bayes CL and Vd of this model as features; only the popPK model is expressible as an nlmixr2 model."
  reference <- "Ma P, Ma H, Liu R, Wen H, Li H, Huang Y, Li Y, Xiong L, Xie L, Wang Q. Prediction of vancomycin plasma concentration in elderly patients based on multi-algorithm mining combined with population pharmacokinetics. Sci Rep. 2024;14(1):27165. doi:10.1038/s41598-024-78558-1"
  vignette <- "Ma_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ma 2024 Results, 'Population pharmacokinetic model': 'CLcr was identified as the only significant",
        "covariate for CL' and 'CLcr was calculated according to the Cockroft-Gault Equations'.",
        "The Results CL equation and Table 2 both normalize by 65.24 mL/min; the paper does not state what",
        "cohort statistic 65.24 is, and Table 1 tabulates eGFR (median 86.85 mL/min/1.73 m^2), not CLcr, so",
        "the CLcr distribution of the cohort is not recoverable from the paper. Stored under canonical CRCL",
        "with raw mL/min, following the inst/references/covariate-columns.md CRCL entry (which accepts a raw",
        "Cockcroft-Gault creatinine clearance when the source paper does not BSA-normalize) and the",
        "precedent in Delattre_2010_amikacin.R and Alqahtani_2018_vancomycin.R. Do NOT feed a",
        "BSA-normalized eGFR into this model in place of CRCL: the two scales differ systematically in an",
        "elderly cohort and the covariate enters as a power term.",
        "NOTE (see the vignette Errata): the published article prints 'CLcr was calculated according to the",
        "Cockroft-Gault Equations24:' followed by NO equation -- the display equation is absent from the",
        "typeset article. Confirmed three ways: rendered at 300 dpi, page 5 shows the colon followed by",
        "white space and the PDF carries no image object on that page; the publisher's per-equation asset",
        "set served by the EuropePMC supplementaryFiles endpoint contains exactly two equation graphics",
        "(Article_Equa = the CL equation, Article_Equb = the V equation) and none for Cockcroft-Gault; and",
        "the Supplementary Information (MOESM1) is figure captions S1-S4 only. The standard Cockcroft-Gault",
        "form cited as reference 24 must therefore be assumed when generating CLcr for this model."
      ),
      source_name        = "CLcr"
    )
  )

  # Screened in the Ma 2024 stepwise covariate analysis but NOT retained in the
  # final model. Methods, 'Population pharmacokinetic analysis': "Covariates
  # included all features listed in Table 1, excluding vancomycin
  # concentration, vancomycin administration, and pop PK parameters"; Results,
  # 'Population pharmacokinetic model': "CLcr was identified as the only
  # significant covariate for CL". No point estimate is published for any of
  # these, so they are documentation only. Table 1 also screened covariates
  # that have no canonical column in inst/references/covariate-columns.md and
  # are therefore recorded only in population$notes: SOFA score, dialysis
  # (CRRT / PD / HD), uric acid, respiratory failure type, co-medication count,
  # diabetes, hypertension, hyperlipidemia, and procalcitonin.
  covariatesDataExcluded <- list(
    AGE   = list(description = "Age", units = "years", type = "continuous",
                 notes = "Ma 2024 Table 1: median 70 (IQR 66-75) years, training group. Screened; not retained."),
    SEXF  = list(description = "Female sex indicator", units = "binary (1 = female)", type = "categorical",
                 notes = "Ma 2024 Table 1: 58.50% female in the training group (source column 'Gender'). Screened; not retained."),
    BMI   = list(description = "Body mass index (recorded categorically)", units = "kg/m^2", type = "categorical",
                 notes = "Ma 2024 Table 1: <18.5 in 10.78%, 18.5-23.9 in 66.99%, >=24 in 22.22% (training). 26.37% missing overall, median-imputed. Screened; not retained."),
    CREAT = list(description = "Serum creatinine", units = "umol/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 71 (IQR 52.8-105.55) umol/L, training group. Screened; not retained in the popPK model, but retained as a feature of the machine-learning ensemble."),
    CYSC  = list(description = "Cystatin C", units = "mg/L", type = "continuous",
                 notes = "Ma 2024: 39.16% missing overall, median-imputed. Screened; not retained in the popPK model, but retained as a feature of the machine-learning ensemble."),
    ALB   = list(description = "Serum albumin", units = "g/L", type = "continuous",
                 notes = "Ma 2024 Table 1: mean 31.74 (SD 4.23) g/L, training group. Screened; not retained."),
    TPRO  = list(description = "Total protein", units = "g/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 61.1 (IQR 56.08-66.5) g/L, training group. Screened; not retained in the popPK model, but retained as a feature of the machine-learning ensemble."),
    ALT   = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 17.9 (IQR 10.98-34.1) U/L, training group. Screened; not retained in the popPK model, but retained as a feature of the machine-learning ensemble."),
    AST   = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 25.7 (IQR 17.58-39.43) U/L, training group. Screened; not retained."),
    GGT   = list(description = "Gamma-glutamyltransferase", units = "U/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 39 (IQR 22.7-93.3) U/L, training group. Screened; not retained."),
    ALP   = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 99 (IQR 74-151) U/L, training group. Screened; not retained."),
    TBILI = list(description = "Total bilirubin", units = "umol/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 12.42 (IQR 9.4-20) umol/L, training group. Screened; not retained."),
    WBC   = list(description = "White blood cell count", units = "10^9/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 7.62 (IQR 5.57-11.02) x10^9/L, training group. Screened; not retained."),
    HGB   = list(description = "Hemoglobin", units = "g/L", type = "continuous",
                 notes = "Ma 2024 Table 1: median 88 (IQR 78-104.25) g/L, training group. Screened; not retained in the popPK model, but retained as a feature of the machine-learning ensemble."),
    NEUT  = list(description = "Neutrophil percentage", units = "% of white blood cells", type = "continuous",
                 notes = "Ma 2024 Table 1 source column 'NEU%': median 72.6 (IQR 61.48-83.48)%, training group. Screened; not retained.")
  )

  population <- list(
    species          = "human",
    n_subjects       = 245L,
    n_studies        = 1L,
    age_range        = ">=65 years (inclusion criterion)",
    age_median       = "70 years (IQR 66-75), training group",
    weight_range     = "Not reported; body mass index tabulated categorically only",
    weight_median    = "Not reported",
    sex_female_pct   = 58.5,
    race_ethnicity   = "Not reported (single-center study in Chongqing, China)",
    disease_state    = "Elderly inpatients (age >=65 years) with suspected or documented Gram-positive bacterial infection, treated with vancomycin for >=2 days and undergoing therapeutic drug monitoring. Comorbidities tabulated in Table 1 include diabetes (59.15%), hypertension (46.73%), respiratory failure (25.82%), and hyperlipidemia (2.6%); 21.20% received dialysis (CRRT, peritoneal dialysis, or hemodialysis) and 43.14% had a SOFA score >=2.",
    dose_range       = "Total daily dose median 1500 mg (IQR 800-2000) in the training group and 1000 mg (IQR 1000-2000) in the testing group; single dose median 13.33 mg/kg (IQR 8.89-16.71). Dosing interval 12 h in 66.10%, 24 h in 16.99%, other in 16.99%. Intravenous infusion in 100 mL (37.25%) or 250 mL (55.88%) diluent; infusion duration is not reported. Median total treatment duration 5 days (IQR 3-10).",
    regions          = "China (single center, Southwest Hospital / First Affiliated Hospital of Army Medical University, Chongqing)",
    renal_function   = "eGFR median 86.85 mL/min/1.73 m^2 (IQR 53.97-103.87), training group. The model's CLcr covariate is a raw Cockcroft-Gault creatinine clearance normalized to 65.24 mL/min; the CLcr distribution itself is not tabulated.",
    n_concentrations = 383L,
    notes            = "Single-center retrospective study, November 2013 to July 2022 (Ma 2024 Methods, 'Patients and data'; Table 1). 383 TDM measurements from 245 elderly patients, split 8:2 into a training group (n = 306 measurements) and a testing group (n = 77 measurements); Table 1 percentages are per measurement, not per patient. Samples were drawn within 30 min before the next morning dose after at least two days of continuous administration, so essentially every observation is a steady-state trough (observed concentration median 15.1 mg/L, IQR 10.7-20.63, in the training group and 16.9 mg/L, IQR 11.95-21.9, in the testing group; time since last dose median 10.89 h). Assay: enzyme-multiplied immunoassay technique (EMIT) on a Viva-ProE system (Syva, USA). Model fit in NONMEM 7.5.1 with FOCE-I; evaluated by bootstrap (981 of 1000 datasets converged) and VPC (Figure 2). Covariates screened in the stepwise analysis but not retained are listed in covariatesDataExcluded, plus the following which have no canonical covariate column: SOFA score, dialysis (CRRT / peritoneal dialysis / hemodialysis), uric acid, respiratory failure type, co-medication count, diabetes, hypertension, hyperlipidemia, and procalcitonin. Three further clinical validation cohorts (27, 40, and 25 patients) were used to evaluate the machine-learning ensemble, not to refit the popPK model; the paper dates their enrolment inconsistently, as 'August 2022 to May 2023' in Methods 'Patients and data' and as 'August 2022 to September 2024' in Methods 'Modeling and validation'. Neither window affects the popPK model, which was fit to the 383 measurements above."
  )

  ini({
    # Structural parameters (Ma 2024 Table 2 and the Results 'Population
    # pharmacokinetic model' equations, which read in full:
    #   CL (L/h) = 3.02 * (CLcr/65.24)^0.856 * exp(etaCL)
    #   V  (L)   = 83.3 * exp(etaV)
    # A one-compartment model with first-order elimination; no absorption
    # phase because vancomycin is given as an intravenous infusion.
    lcl <- log(3.02); label("Clearance at CRCL = 65.24 mL/min (L/h)")  # Ma 2024 Table 2 theta_1 = 3.02 (RSE 2.7%; bootstrap median 3.00, 95% CI 2.85-3.17)
    lvc <- log(83.3); label("Volume of distribution (L)")              # Ma 2024 Table 2 theta_3 = 83.3 (RSE 8.5%; bootstrap median 83.7, 95% CI 70.1-101)

    # Covariate effect (Ma 2024 Table 2 row "CL = theta_1 * (CLcr/65.24)^theta_2").
    # theta_2 is the "typical value of creatinine clearance-dependent fraction
    # of clearance" per the Table 2 legend.
    e_crcl_cl <- 0.856; label("Power exponent on (CRCL/65.24) for CL (unitless)")  # Ma 2024 Table 2 theta_2 = 0.856 (RSE 4.1%; bootstrap median 0.857, 95% CI 0.792-0.936)

    # Inter-individual variability. Methods, 'Population pharmacokinetic
    # analysis': "To estimate the interindividual variability of the
    # pharmacokinetic parameters, an exponential model was applied", and the
    # Results equations write the random effects as exp(etaCL) / exp(etaV).
    # Table 2 reports the two IIV rows under the header "Inter-individual
    # variability (%)", i.e. as coefficients of variation, so for a log-normal
    # eta the variance is omega^2 = log(1 + CV^2):
    #   CL: log(1 + 0.308^2) = log(1.094864) = 0.090630
    #   V:  log(1 + 0.488^2) = log(1.238144) = 0.213613
    # The alternative reading -- that the printed percentages ARE omega x 100,
    # giving 0.094864 and 0.238144 -- shifts the CL SD by 1.9% and the V SD by
    # 5.3% relative, and the paper gives no signal that discriminates the two
    # (the proportional residual row is numerically identical on both scales).
    # The CV reading is used here because it matches the exponential-IIV
    # declaration in Methods and the convention applied across the other
    # vancomycin models in this package. See the vignette Errata.
    etalcl ~ 0.090630  # Ma 2024 Table 2: IIV_CL = 30.8% (RSE 8.8%; bootstrap median 30.4, 95% CI 23.3-35.5); log(1 + 0.308^2)
    etalvc ~ 0.213613  # Ma 2024 Table 2: IIV_V  = 48.8% (RSE 17.8%; bootstrap median 47.5, 95% CI 19.3-66.1); log(1 + 0.488^2)

    # Combined additive plus proportional residual error. Results, 'Population
    # pharmacokinetic model': "The combined model was selected to evaluate the
    # residual variability."
    propSd <- 0.208; label("Proportional residual error (fraction)")  # Ma 2024 Table 2: Prop_error = 20.8% (RSE 14.7%; bootstrap median 20.1, 95% CI 13.5-28.2)
    addSd  <- 1.92;  label("Additive residual error (mg/L)")          # Ma 2024 Table 2: Add_error = 1.92 mg/L (RSE 39.1%; bootstrap median 1.89, 95% CI 0.826-3.31)
  })
  model({
    # Individual PK parameters. CL scales by (CRCL/65.24)^0.856; V has no
    # retained covariate.
    cl <- exp(lcl + etalcl) * (CRCL / 65.24)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
