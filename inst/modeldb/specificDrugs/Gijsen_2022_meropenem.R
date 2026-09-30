Gijsen_2022_meropenem <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "intermittently infused meropenem (30-minute infusions) in critically ill",
    "adults with severe sepsis or septic shock and preserved or increased",
    "renal function (Gijsen 2022; 58 patients, 345 plasma concentrations over",
    "70 dosing intervals, eGFR CKD-EPI >= 70 mL/min/1.73 m^2, no renal",
    "replacement therapy). Clearance (13.7 L/h at the reference) scales as a",
    "power function of Cockcroft-Gault creatinine clearance normalised to",
    "111.7 mL/min (exponent 0.637); no other covariate was retained.",
    "Correlated interindividual variability on clearance and central volume,",
    "independent variability on peripheral volume, and combined proportional",
    "plus additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Gijsen M, Elkayal O, Annaert P, Van Daele R, Meersseman P, Debaveye Y,",
    "Wauters J, Dreesen E, Spriet I (2022). Meropenem Target Attainment and",
    "Population Pharmacokinetics in Critically Ill Septic Patients with",
    "Preserved or Increased Renal Function. Infection and Drug Resistance",
    "15:53-62. doi:10.2147/IDR.S343264.",
    "Final parameter estimates from the NONMEM control stream in",
    "Supplementary File S2, cross-checked against Table 2.",
    sep = " "
  )
  vignette <- "Gijsen_2022_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source name CG in the File S2 $INPUT block, eCrClCG in the text.",
        "Raw mL/min, NOT BSA-normalized: Equation 2 prints the reference as",
        "'111.7 mL/min', matching the raw-mL/min variant of the canonical",
        "CRCL column (Delattre_2010_amikacin.R / Jin_2026_colistinSulfate.R",
        "precedents). Enters clearance as the power function",
        "(CRCL / 111.7)^0.637 per the File S2 $PK line",
        "'CLCG = ((CG/111.7) **THETA (5))'. The paper does not state what",
        "statistic 111.7 mL/min is; it is not tabulated in Table 1, which",
        "reports only measured 24-hour urinary creatinine clearance (median",
        "84 and 109 mL/min/1.73 m^2 on early and late sampling days) and",
        "CKD-EPI eGFR. It is treated as the typical-patient value. The",
        "Cockcroft-Gault body-weight convention (total, ideal or adjusted",
        "body weight) is not stated. Results: 'eCrClCG was the only",
        "covariate withheld in our final model, explaining 18% of the IIV on",
        "CL. The eCrClCG performed significantly better than any other",
        "covariate.'"
      ),
      source_name = "CG"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 'Population PK Modelling') and not retained. Cohort median 63 [IQR 55; 68] years (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened as a covariate and, separately, as allometric scaling",
        "(exponents 0.75 on CL and 1 on volumes); neither was retained",
        "('Allometric scaling was not retained since it did not improve the",
        "base model'). Cohort median 70 [IQR 60; 79] kg (Table 1)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained. Cohort distribution not reported."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained. Cohort distribution not reported."
    ),
    ABW = list(
      description = "Adjusted body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained. Cohort distribution not reported."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained. Cohort distribution not reported."
    ),
    RENAL_ARC = list(
      description = "Augmented renal clearance indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "The screened covariate was the PROBABILITY of augmented renal",
        "clearance from the ARC predictor (Methods reference 30), not a",
        "binary indicator; it was not retained. ARC (measured 24-hour urinary",
        "creatinine clearance >= 130 mL/min/1.73 m^2) was present on 16 of 66",
        "dosing intervals with a measurement (24.2%)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    n_observations = 345L,
    n_dosing_intervals = 70L,
    age_median = "63 years (IQR 55-68)",
    weight_median = "70 kg (IQR 60-79)",
    sex_female_pct = 100 * 18 / 58,
    disease_state = paste(
      "Critically ill adults in the intensive care unit with severe sepsis",
      "(45%) or septic shock (55%), mainly respiratory infection focus (72%).",
      "APACHE II 20 [16; 25] and SOFA 9 [7; 10] on ICU admission. Excluded:",
      "pregnancy, do-not-resuscitate order, extracorporeal membrane",
      "oxygenation, renal replacement therapy, and eGFR CKD-EPI < 70",
      "mL/min/1.73 m^2 on the sampling day."
    ),
    renal_function = paste(
      "Preserved or increased. Measured 24-hour urinary creatinine clearance",
      "84 [64; 119] and 109 [75; 136] mL/min/1.73 m^2 on early and late",
      "sampling days (overall 88 [69; 128]); CKD-EPI eGFR 104 [91; 119] and",
      "112 [101; 119] mL/min/1.73 m^2. Augmented renal clearance on 24.2% of",
      "dosing intervals. The Cockcroft-Gault covariate distribution is not",
      "tabulated; the model reference is 111.7 mL/min."
    ),
    dose_range = paste(
      "Meropenem 1000 mg q8h on 64% of sampling days, 2000 mg q8h on most of",
      "the rest, 500 mg q8h on two sampling days; all as 30-minute",
      "intravenous infusions. Daily dose 3000 [3000; 6000] mg."
    ),
    regions = "Belgium (University Hospitals Leuven, single centre).",
    notes = paste(
      "Prospective observational cohort, October 2013 to October 2017",
      "(NCT03560557). Rich sampling around one dosing interval on an early",
      "(day 2 +/- 1) and/or late (day 5 +/- 1) day of therapy at predose, 30,",
      "120 and 240 minutes after infusion start and 15 minutes before the next",
      "dose. 83% of patients sampled on an early day, 12 patients on both.",
      "Total concentrations were treated as unbound (2% protein binding).",
      "UPLC-MS/MS assay, LLOQ 0.09 mg/L; no sample fell below it. Samples",
      "stored at -20 C in year one were corrected with a degradation model.",
      "NONMEM 7.4, FOCE with interaction, ADVAN13. Inter-occasion variability",
      "could not be identified. Demographics from Table 1."
    )
  )

  ini({
    # Structural parameters -- Supplementary File S2 control stream $THETA
    # (final estimates; they agree with Table 2 to its printed precision).
    lcl <- log(13.7); label("Clearance at CRCL = 111.7 mL/min (L/h)")  # File S2 THETA(1) = 13.7; Table 2 CL = 13.7 L/h (RSE 7.4%)
    lvc <- log(25.5); label("Central volume of distribution (L)")      # File S2 THETA(2) = 25.5; Table 2 Vc = 25.5 L (RSE 8.5%)
    lq <- log(8.13); label("Intercompartmental clearance (L/h)")       # File S2 THETA(3) = 8.13; Table 2 Q = 8.1 L/h (RSE 32.6%)
    lvp <- log(12.4); label("Peripheral volume of distribution (L)")   # File S2 THETA(4) = 12.4; Table 2 Vp = 12.4 L (RSE 15.8%)

    # Power exponent of Cockcroft-Gault CrCl on clearance. File S2 THETA(5) =
    # 0.637 and Table 2 'eCrCl CG on CL' = 0.64 (bootstrap median 0.64) agree;
    # Equation 2 in the Results prints 0.725, a typo (see vignette Errata).
    e_crcl_cl <- 0.637; label("Power exponent of CRCL on clearance (unitless)")  # File S2 THETA(5) = 0.637; Table 2 = 0.64 (RSE 20.9%)

    # Interindividual variability. File S2 $OMEGA BLOCK(2) on CL and V1 plus a
    # separate $OMEGA on V2, all exponential (Equation 1: CL_i = TVCL x
    # exp(eta_i)). The values are variances; Table 2 reports them as
    # CV% = sqrt(exp(omega^2) - 1) x 100: sqrt(exp(0.288)-1) = 57.8%,
    # sqrt(exp(0.281)-1) = 57.0%, sqrt(exp(0.441)-1) = 74.4%. The BLOCK(2)
    # off-diagonal 0.185 is a COVARIANCE (correlation 0.185 / sqrt(0.288 x
    # 0.281) = 0.65); Table 2 labels that same number 'Correlation between
    # CL & Vc'. No IIV on Q.
    etalcl + etalvc ~ c(0.288, 0.185, 0.281)  # File S2 $OMEGA BLOCK(2): 0.288 (IIV CL), 0.185 (cov), 0.281 (IIV V1)
    etalvp ~ 0.441                            # File S2 $OMEGA: 0.441 (IIV V2); Table 2 74.4% CV

    # Residual error. File S2 $ERROR: Y = IPRED * (1 + EPS(1)) + EPS(2) with
    # a diagonal $SIGMA of 0.147 (EPS(1), commented 'prop err') and 0.0925
    # (EPS(2), commented 'add err'). These are variances, so the SDs are
    # sqrt(0.147) = 0.383 (proportional) and sqrt(0.0925) = 0.304 mg/L
    # (additive). Table 2 swaps the labels: its 'Proportional residual
    # variability 31.1 %CV' is sqrt(exp(0.0925) - 1) = 31.1% computed from the
    # ADDITIVE variance, and its 'Additive residual variability 0.147 mg/L' is
    # the raw PROPORTIONAL variance. The as-run $ERROR code is encoded here.
    propSd <- 0.3834; label("Proportional residual error (fraction)")  # File S2 $SIGMA(1,1) = 0.147 (variance, EPS(1) proportional) -> sqrt = 0.3834
    addSd <- 0.3041; label("Additive residual error (mg/L)")           # File S2 $SIGMA(2,2) = 0.0925 (variance, EPS(2) additive) -> sqrt = 0.3041
  })

  model({
    # 1. Individual PK parameters (File S2 $PK). Clearance scales with
    # Cockcroft-Gault CrCl by a power function normalised to 111.7 mL/min
    # (CLCG = (CG/111.7)**THETA(5); TVCL = CLCOV * THETA(1)).
    cl <- exp(lcl + etalcl) * (CRCL / 111.7)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    # 2. Micro-constants (File S2: KE = CL/V1, K12 = Q/V1, K21 = Q/V2).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with linear elimination (File S2 $DES).
    # Meropenem was given as 30-minute intravenous infusions into the central
    # compartment (DEFDOSE); there is no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. Plasma meropenem in mg/L (S1 = V1); total
    # concentrations were treated as unbound (2% protein binding). Combined
    # proportional plus additive residual error with variances summed (File
    # S2: W = SQRT(SIGMA(1,1) * IPRED**2 + SIGMA(2,2))).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
