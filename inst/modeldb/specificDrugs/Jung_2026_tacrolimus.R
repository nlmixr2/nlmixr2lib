Jung_2026_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption (ka fixed) for oral tacrolimus in pediatric liver transplant recipients younger than 36 months during the immediate postoperative period (postoperative days 2-11), fitted to routine whole-blood trough monitoring data. Apparent clearance CL/F falls with total bilirubin and aspartate aminotransferase (negative power effects); apparent volume V/F rises with age (fixed power exponent). Correlated inter-individual variability on CL/F and V/F, plus inter-occasion variability on CL/F with a separately estimated variance for each of six postoperative occasions; combined additive and proportional residual error."
  reference <- paste(
    "Jung SM, Kim M, Lee SH, Choi S, Bae KS, Lim HS (2026).",
    "Population pharmacokinetic modeling of tacrolimus using real-world data for",
    "precision dosing in pediatric liver transplant recipients younger than 36",
    "months during the immediate postoperative period.",
    "Front Pharmacol 17:1896778. doi:10.3389/fphar.2026.1896778.",
    "Parameter values are Tables 4 and 5 of the paper; the covariate equations",
    "are Equations 9 and 10.",
    sep = " "
  )
  vignette <- "Jung_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    TBILI = list(
      description = "Total serum bilirubin (time-varying)",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Jung 2026 reports and models total bilirubin in mg/dL with reference",
        "2 mg/dL: TV(CL/F) = 4.85 * (TBIL/2)^-0.353 * (AST/60)^-0.106",
        "(Equation 9, Table 4). The canonical column is SI umol/L, so model()",
        "converts it back with TBILI / 17.1 before applying the published",
        "reference. Time-varying; laboratory values not measured on the day of",
        "a tacrolimus sample were carried forward from the most recent prior",
        "measurement (Methods). Cohort subject-level median 2.30 mg/dL (range",
        "0.80-10.60); individual records 0.4-29.0 mg/dL (Table 1). Occasion",
        "medians fell from 4.60 mg/dL on POD 2 to 1.20 mg/dL on POD 9-11",
        "(Table 7).",
        sep = " "
      ),
      source_name = "TBIL"
    ),
    AST = list(
      description = "Aspartate aminotransferase (time-varying)",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 60 U/L, negative power exponent 0.106 on CL/F (Equation 9,",
        "Table 4). Time-varying, last-observation-carried-forward. Cohort",
        "subject-level median 68 U/L (range 26-461); individual records",
        "17-8,781 U/L (Table 1). Occasion medians fell from 188 U/L on POD 2 to",
        "45 U/L on POD 9-11 (Table 7). AST was retained in preference to ALT",
        "because the two are collinear (Results, Covariate screening).",
        sep = " "
      ),
      source_name = "AST"
    ),
    AGE = list(
      description = "Recipient age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Jung 2026 models age in months with reference 12 months:",
        "TV(V/F) = 124 * (AGE/12)^0.16 (Equation 10; exponent fixed at its",
        "forward-selection estimate, Table 4). The canonical AGE column is in",
        "years, so model() multiplies by 12 before applying the published",
        "reference; a 1-year-old is the reference individual. Cohort median",
        "11.70 months (range 3.30-35.50) (Table 1).",
        sep = " "
      ),
      source_name = "AGE"
    ),
    OCC = list(
      description = "Pharmacokinetic occasion (OCCPK) index for inter-occasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Jung 2026 OCCPK, defined by postoperative day (POD): OCC 1 = POD 2,",
        "2 = POD 3, 3 = POD 4, 4 = POD 5-6, 5 = POD 7-8, 6 = POD 9-11",
        "(Methods, Fixed effect and random effect; Table 2). Each occasion has",
        "its own estimated IOV variance (Table 5), unlike the shared",
        "$OMEGA BLOCK(1) SAME structure of most IOV models. Tacrolimus was",
        "started on POD 1, before the first occasion; the paper does not say",
        "which occasion its dosing records carried, and the dosing records of",
        "POD 1 are assigned OCC = 1 here (the CL/F that governs accumulation",
        "into the POD 2 trough). OCC = 0 or any value outside 1-6 zeros every",
        "indicator and gives CL/F without IOV; the model was not fitted beyond",
        "POD 11.",
        sep = " "
      ),
      source_name = "OCCPK"
    )
  )

  # Covariates screened by Jung 2026 (Supplementary Tables S1 and S2) but not
  # retained in the final model. Graft weight, age at operation and ABO blood
  # type (as a single categorical) were also screened and not retained; they
  # have no single canonical column and are recorded here in prose only.
  # Graft weight: CL/F exponent 0.249 (dOFV -1.66), V/F -0.137 (dOFV -0.58).
  # Age at operation: CL/F 0.021 (dOFV -0.035), V/F 0.108 (dOFV -1.18).
  # ABO type: -0.001 on both (dOFV -0.001).
  covariatesDataExcluded <- list(
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Second-largest univariate CL/F effect (exponent -0.290, dOFV -62.88, Table S1) and a further dOFV -10.13 on top of TBIL (Table S2), but dropped in favour of AST because the two transaminases are collinear.",
      source_name = "ALT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F exponent 0.453 on top of TBIL (dOFV -11.32, Table S2), not retained because exogenous albumin administration confounds it in the early postoperative period.",
      source_name = "ALB"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F exponent -0.197 on top of TBIL (dOFV -5.24, Table S2); not retained on statistical and pharmacological grounds (tacrolimus is hepatically cleared and itself raises creatinine).",
      source_name = "SCR"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F exponent -0.221 (dOFV -1.30, Table S1); V/F exponent -0.149 (dOFV -1.90, Table S2), a biologically implausible negative exponent, so excluded.",
      source_name = "WT"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F exponent -0.108 (dOFV -0.42, Table S1). Not retained.",
      source_name = "Hematocrit"
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F exponent 0.047 (dOFV -0.17); V/F 0.017 (dOFV -0.04) (Table S1). Not retained.",
      source_name = "Hemoglobin"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Categorical shift 0.022 on CL/F and V/F (dOFV -0.04, Table S1). Not retained.",
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_observations = 731L,
    n_studies = 1L,
    age_range = "3.3-35.5 months (median 11.7 months); age at surgery 3.2-35.4 months",
    age_median = "11.7 months",
    weight_range = "4.3-14.4 kg",
    weight_median = "8.95 kg",
    sex_female_pct = 58.7,
    race_ethnicity = "Korean (single centre, Seoul); not tabulated",
    disease_state = "Pediatric recipients of a primary living- or deceased-donor liver transplant for biliary atresia, younger than 36 months at transplantation, observed during the immediate postoperative period (POD 2-11).",
    dose_range = "Oral tacrolimus twice daily from postoperative day 1, initial dose 0.0375 mg/kg per dose per institutional protocol, then adjusted by routine trough monitoring (target trough 10 ng/mL).",
    regions = "South Korea (Asan Medical Center Children's Hospital)",
    notes = "Retrospective single-centre cohort, January 2010 - May 2023. 731 whole-blood tacrolimus concentrations (predominantly morning troughs, ACMIA immunoassay, LLOQ 1.0 ng/mL) from 92 subjects, collected through POD 11 and grouped into six occasions (Table 2). Table 1: graft weight 287.5 g (120-915), hematocrit 27.45% (21.4-34.1), albumin 3.15 g/dL (2.2-4.75), creatinine 0.17 mg/dL (0.10-1.70), all subject-level medians (range). Estimated in NONMEM 7.6 by FOCE-I."
  )

  ini({
    # Structural parameters: Jung 2026 Table 4.
    lka <- fixed(log(4.5))
    label("Absorption rate constant ka (1/h)") # Table 4, KA = 4.5 h-1 (fixed); Methods: taken from Wang et al. 2019 pediatric tacrolimus
    lcl <- log(4.85)
    label("Apparent clearance CL/F at TBIL 2 mg/dL and AST 60 U/L (L/h)") # Table 4, TV(CL/F) = 4.85 L/h (RSE 6.0%); Equation 9
    # Table 4 prints 123.9 and Equation 10 rounds it to 124; the typical V/F
    # values of Table 8 (e.g. 139.0 L at 24.6 months) reproduce only with 123.9.
    lvc <- log(123.9)
    label("Apparent volume of distribution V/F at age 12 months (L)") # Table 4, TV(V/F) = 123.9 L (RSE 9.9%); Equation 10

    # Covariate effects. Table 4 calls TBIL and AST 'negative power' effects
    # with the printed exponents positive; Equation 9 carries the signs, so the
    # stored values are the signed exponents.
    e_tbili_cl <- -0.353
    label("Power exponent of (TBIL / 2 mg/dL) on CL/F (unitless)") # Table 4, theta TBIL = 0.353 (RSE 14.0%) as a negative power; Equation 9
    e_ast_cl <- -0.106
    label("Power exponent of (AST / 60 U/L) on CL/F (unitless)") # Table 4, theta AST = 0.106 (RSE 46.7%) as a negative power; Equation 9
    e_age_vc <- fixed(0.16)
    label("Power exponent of (age / 12 months) on V/F (unitless)") # Table 4, theta AGE = 0.16 (fixed, from forward selection 0.156 in Table S2); Equation 10

    # Inter-individual variability (Table 5): variances, full block, r = 0.7.
    # Covariance = 0.7 * sqrt(0.319 * 0.258) = 0.2008.
    etalcl + etalvc ~ c(0.319, 0.2008, 0.258) # Table 5, IIV CL/F 0.319 (CV 61.3%), V/F 0.258 (CV 54.2%), correlation r 0.7

    # Inter-occasion variability on CL/F (Table 5): one estimated variance per
    # occasion, none shared.
    etaiov_cl_1 ~ 0.929 # Table 5, IOV CL/F OCCPK 1 (POD 2) 0.929 (CV 123.8%)
    etaiov_cl_2 ~ 0.114 # Table 5, IOV CL/F OCCPK 2 (POD 3) 0.114 (CV 34.8%)
    etaiov_cl_3 ~ 0.034 # Table 5, IOV CL/F OCCPK 3 (POD 4) 0.034 (CV 18.7%)
    etaiov_cl_4 ~ 0.011 # Table 5, IOV CL/F OCCPK 4 (POD 5-6) 0.011 (CV 10.7%)
    etaiov_cl_5 ~ 0.044 # Table 5, IOV CL/F OCCPK 5 (POD 7-8) 0.044 (CV 21.2%)
    etaiov_cl_6 ~ 0.14 # Table 5, IOV CL/F OCCPK 6 (POD 9-11) 0.14 (CV 38.8%)

    # Residual error (Table 4; Equation 3 gives W = sqrt(add^2 + prop^2 * IPRED^2),
    # so both are standard deviations).
    addSd <- 1.37
    label("Additive residual error (ng/mL)") # Table 4, sigma add = 1.37 ng/mL (RSE 19.0%)
    propSd <- 0.206
    label("Proportional residual error (fraction)") # Table 4, sigma prop = 0.206 (RSE 7.8%)
  })

  model({
    # Occasion indicators for the per-occasion IOV etas on CL/F. OCC outside
    # 1-6 zeros every indicator.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # Covariates on the paper's scales: bilirubin in mg/dL, age in months.
    tbili_mgdl <- TBILI / 17.1
    age_months <- AGE * 12

    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl) * (tbili_mgdl / 2)^e_tbili_cl * (AST / 60)^e_ast_cl
    vc <- exp(lvc + etalvc) * (age_months / 12)^e_age_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L; x1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
