Fu_2022_teicoplanin <- function() {
  description <- paste(
    "One-compartment IV-infusion population PK model for UNBOUND (ultrafiltrate)",
    "teicoplanin in 72 Chinese adult inpatients (Fu 2022). The central state",
    "holds unbound drug only, so V is the apparent volume of distribution of",
    "unbound teicoplanin (811 L, about ten times the total-drug volume) and CL",
    "is unbound clearance. CKD-EPI estimated glomerular filtration rate enters",
    "CL as a power term normalised to 84 and serum albumin enters V as a power",
    "term normalised to 32 g/L. Exponential between-subject variability on CL",
    "and V; combined proportional + additive residual error. The paper uses",
    "the model in Monte Carlo simulations to recommend eGFR- and",
    "albumin-stratified loading (q12h x 3) and maintenance (q24h) regimens",
    "against unbound trough targets of 0.75 and 1.13 mg/L."
  )
  reference <- paste(
    "Fu W-Q, Tian T-T, Zhang M-X, Song H-T, Zhang L-L. Population",
    "pharmacokinetics and dosing optimization of unbound teicoplanin in Chinese",
    "adult patients. Front Pharmacol. 2022;13:1045895.",
    "doi:10.3389/fphar.2022.1045895"
  )
  vignette <- "Fu_2022_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "teicoplanin, unbound", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the creatinine-based CKD-EPI equation (power effect on unbound clearance)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column eGFR. Fu 2022 Methods Eqs 1-4 print the 2009",
        "creatinine-based CKD-EPI equation (144 or 141 x (Scr/0.7 or 0.9)^a x",
        "0.993^Age, with Scr in mg/dL), which returns a BSA-normalised value in",
        "mL/min/1.73 m^2. The paper labels the column 'ml/min' throughout",
        "(Table 1, Table 3, Methods 'Simulation and dosing optimization'), but no",
        "de-normalisation by body surface area is described and height was not",
        "collected, so the value supplied is the equation output on its native",
        "BSA-normalised scale, which is the canonical CRCL scale. Enters CL as",
        "(CRCL / 84)^0.476 (Results Eq 10). The normalising constant 84 is the",
        "covariate median per the Methods covariate equation (Eq 7, cov_median);",
        "Table 1 reports the mean 80.8 +/- 27.1 (range 17.4-134.0). The paper's",
        "simulations span 20-130."
      ),
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin (power effect on the unbound volume of distribution)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column ALB, g/L (Table 1: 31.0 +/- 5.1 g/L, range 18.3-46.0).",
        "Enters V as (ALB / 32)^1.6 (Results Eq 11), with 32 g/L the covariate",
        "median (Methods Eq 7). The POSITIVE exponent is the expected direction",
        "for an UNBOUND volume: teicoplanin is 90-95% albumin-bound, so lower",
        "albumin raises the unbound fraction and the unbound concentration for a",
        "given total amount, which the model expresses as a smaller apparent",
        "unbound volume (Discussion). The paper's simulations span 15-40 g/L."
      ),
      source_name = "ALB"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "69 +/- 20 years (range 18-99), Fu 2022 Table 1. Screened in the stepwise covariate search and not retained (Results, 'Population pharmacokinetics modeling')."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "26 of 72 patients female (Fu 2022 Table 1). Screened and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "61.0 +/- 9.8 kg (range 40-85), Fu 2022 Table 1. Screened and not retained; the Discussion attributes this to a narrow weight distribution (over 80% between 50 and 70 kg), inaccurate weights in long-term bedridden patients, and the small sample."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "86.2 +/- 48.3 umol/L (Fu 2022 Table 1). Screened and not retained; renal function enters through the CKD-EPI eGFR instead."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "9.5 +/- 7.4 mmol/L (Fu 2022 Table 1). Screened and not retained."
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      notes = "1.4 +/- 0.8 mg/L (Fu 2022 Table 1). Screened and not retained."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "10.6 +/- 5.3 x 10^9/L (Fu 2022 Table 1). Screened and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 72L,
    n_studies = 1L,
    n_concentrations = 103L,
    age_range = "18-99 years (mean 69 +/- 20; Fu 2022 Table 1)",
    weight_range = "40-85 kg (mean 61.0 +/- 9.8; Fu 2022 Table 1)",
    sex_female_pct = 36.1,
    race_ethnicity = "Chinese (single centre, 900th Hospital of Joint Logistics Support Force, Fuzhou)",
    disease_state = paste(
      "Adult inpatients (age >= 18 years) treated with intravenous teicoplanin",
      "for Gram-positive infection: pulmonary 58 (80.6%), skin 6, pyemia 4,",
      "bacteraemia 3, abdominal 2, biliary tract 2, liver abscess 1 (Table 1).",
      "Excluded: pregnancy, haemodialysis, disseminated intravascular",
      "coagulation, and continuous renal replacement therapy."
    ),
    dose_range = paste(
      "Teicoplanin by 40-min intravenous infusion; standard regimen 400 mg q12h",
      "for three doses followed by 400 or 200 mg once daily, at physician",
      "discretion; daily dose range 50-1600 mg. Two products (Sanofi-Aventis",
      "44 patients, HISUN 28); product type was screened and not retained."
    ),
    regions = "China (Fuzhou; enrolled January-December 2019)",
    renal_function = "CKD-EPI eGFR 80.8 +/- 27.1 (range 17.4-134.0); serum creatinine 86.2 +/- 48.3 umol/L (Fu 2022 Table 1).",
    albumin = "Serum albumin 31.0 +/- 5.1 g/L (range 18.3-46.0), including hypoalbuminaemic patients below 25 g/L.",
    notes = paste(
      "Prospective single-centre therapeutic-drug-monitoring study. Sparse",
      "sampling, 1-5 samples per patient (mostly troughs, drawn 30 min to 1 h",
      "before the 4th and 6th doses). Unbound concentrations measured in",
      "Centrifree ultrafiltrate (37 C) by UPLC-MS/MS, linear 0.10-8.00 ug/mL;",
      "observed unbound concentrations 1.5 +/- 0.9 ug/mL (range 0.4-4.4).",
      "Estimation: NONMEM 7.2, FOCEI. A two-compartment model did not improve",
      "the fit (dOFV -0.24) and was over-parameterised. Evaluation by",
      "goodness-of-fit plots, a 1000-sample bootstrap (875 successful) and NPDE."
    )
  )

  ini({
    # Fu 2022 Results Eqs 10-11 (final model) and Table 2 'NONMEM Estimates':
    #   CL (L/h) = 11.7 * (eGFR / 84)^0.476 * exp(eta_CL)
    #   V  (L)   = 811  * (ALB / 32)^1.6    * exp(eta_V)
    lcl <- log(11.7); label("Unbound clearance at eGFR = 84 (L/h)")                       # Fu 2022 Table 2: thetaCL = 11.7 L/h (RSE 7.4%; bootstrap median 11.38, 95% CI 9.06-13.4); Eq 10
    lvc <- log(811); label("Unbound volume of distribution at ALB = 32 g/L (L)")           # Fu 2022 Table 2: thetaV = 811 L (RSE 11.1%; bootstrap median 822, 95% CI 616-1025); Eq 11

    e_crcl_cl <- 0.476; label("Power exponent of eGFR on unbound CL (unitless)")           # Fu 2022 Table 2: thetaeGFR = 0.476 (RSE 50%; bootstrap median 0.501, 95% CI 0.229-1.40); Eq 10 exponent
    e_alb_vc <- 1.60; label("Power exponent of serum albumin on unbound V (unitless)")     # Fu 2022 Table 2: thetaALB = 1.60 (RSE 27.4%; bootstrap median 1.55, 95% CI 0.690-2.55); Eq 11 exponent

    # Fu 2022 Table 2 reports IIV as a percentage ('etaCL (%)', 'etaV (%)')
    # under an exponential model (Methods Eq 5), read here as omega x 100
    # (the SD of eta), so omega^2 = (pct / 100)^2. The paper's own Monte Carlo
    # PTA table (Supplementary Table 1, 72 cells) separates this from the
    # log(1 + CV^2) reading: simulated minus published PTA has median
    # +0.05 percentage points under omega x 100 versus +1.35 under the CV
    # reading (see the vignette).
    etalcl ~ 0.148996 # Fu 2022 Table 2: etaCL = 38.6% (RSE 16.5%; shrinkage 35%) -> 0.386^2
    etalvc ~ 0.288369 # Fu 2022 Table 2: etaV = 53.7% (RSE 13.2%; shrinkage 31%) -> 0.537^2

    # Mixed residual error, Methods Eq 6: Y = F * (1 + eps1) + eps2.
    propSd <- 0.183; label("Proportional residual error (fraction)")                      # Fu 2022 Table 2: eps1 = 18.3% (RSE 14.5%)
    addSd <- 0.122; label("Additive residual error (mg/L)")                               # Fu 2022 Table 2: eps2 = 0.122 ug/mL (RSE 41.3%)
  })

  model({
    # CKD-EPI eGFR (canonical CRCL) on CL and serum albumin (g/L) on V, both
    # as power terms normalised to the covariate medians printed in Eqs 10-11.
    cl <- exp(lcl + etalcl) * (CRCL / 84)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (ALB / 32)^e_alb_vc

    kel <- cl / vc

    # 40-min IV infusion into central (Methods, 'Teicoplanin dosing'); the
    # central state is unbound teicoplanin, so Cc is the unbound concentration
    # (mg / L = ug/mL, the assay unit).
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
