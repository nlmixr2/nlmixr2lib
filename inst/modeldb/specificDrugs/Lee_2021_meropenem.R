Lee_2021_meropenem <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for meropenem in",
    "critically ill Korean adults, including patients on extracorporeal",
    "membrane oxygenation (Lee 2021; n = 26, 8 on ECMO, 125 plasma",
    "samples). Total clearance increases linearly with CKD-EPI estimated",
    "glomerular filtration rate centred at 91.57 mL/min/1.73 m^2:",
    "CL = 6.37 * (1 + 0.00925 * (CRCL - 91.57)) L/h. ECMO support was",
    "tested and did not affect any PK parameter. Log-normal IIV on CL, Vc",
    "and Vp (none on Q); residual error is a power model whose standard",
    "deviation is 0.246 * Cc^0.865. The unbound concentration",
    "Cu = fu * Cc (fu = 0.98) drives the fT>MIC targets (40% fT>MIC,",
    "100% fT>MIC, 100% fT>4MIC) of the paper's Monte Carlo",
    "probability-of-target-attainment simulations."
  )
  reference <- paste(
    "Lee DH, Kim HS, Park S, Kim HI, Lee SH, Kim YK. (2021). Population",
    "Pharmacokinetics of Meropenem in Critically Ill Korean Patients and",
    "Effects of Extracorporeal Membrane Oxygenation. Pharmaceutics",
    "13(11):1861. doi:10.3390/pharmaceutics13111861.",
    sep = " "
  )
  vignette <- "Lee_2021_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Lee 2021 Methods 2.4 (plasma meropenem by
  # HPLC-MS/MS) and Results 3.2 (two-compartment model with CL, VC, VP, Q).
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI creatinine equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the Lee 2021 final model (Results",
        "3.2: 'GFR was estimated using the CKD-EPI equation and identified",
        "as a statistically significant covariate of CL'). Enters CL",
        "linearly, CL = theta1 * (1 + theta2 * (CE - 91.57)) (Table 2",
        "structural-model row). The BSA-normalized CKD-EPI form (mL/min/1.73",
        "m^2) is the one used: Table 1 lists it with medians 87.7 (ECMO) and",
        "91.6 (non-ECMO) mL/min/1.73 m^2, bracketing the 91.57 centring",
        "value, and the PTA Figures 3-6 stratify patients by 'eGFR' in",
        "mL/min/1.73 m^2. The BSA-de-normalized 'modified CKD-EPI' variant",
        "(mL/min; medians 82.4 / 77.7) was a separate candidate and was not",
        "the one selected. Among the renal-function candidates the CKD-EPI",
        "equation gave the largest OFV drop on the base model (27.933;",
        "Discussion). Treated as time-fixed per subject."
      ),
      source_name = "CE"
    )
  )

  # Lee 2021 Methods 2.5 screened these covariates by stepwise forward
  # inclusion (p < 0.01) and backward exclusion (p < 0.001); only CKD-EPI
  # eGFR on CL survived. Recorded for provenance; none enters model() and
  # no point estimates are published for them.
  covariatesDataExcluded <- list(
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation support indicator (1 = on ECMO, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving ECMO)",
      notes = paste(
        "The paper's primary question. ECMO did not affect meropenem PK",
        "(Results 3.2; Table 3 compares individual CL, VC, VP and VSS between",
        "the groups, all p > 0.4). ECMO type (veno-arterial n = 7,",
        "veno-venous n = 1) was also screened and not retained."
      ),
      source_name = "ECMO"
    ),
    Q_ECMO = list(
      description = "ECMO circuit blood flow rate",
      units = "L/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as 'ECMO flow rate' (Methods 2.5); not retained. Values are not tabulated.",
      source_name = "ECMO flow rate"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained. Table 1 median (IQR): ECMO 64.0 (56.3-66.5), non-ECMO 72.0 (66.0-80.3) years.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as 'sex' (Methods 2.5); not retained. Table 1: ECMO 4 male / 4 female, non-ECMO 14 male / 4 female.",
      source_name = "Sex"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained. Table 1 median (IQR): ECMO 162 (153-169), non-ECMO 165 (156-170) cm.",
      source_name = "Height"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained, so the model carries no allometric scaling. Table 1 median (IQR): ECMO 63.5 (61.9-66.3), non-ECMO 54.4 (50.5-64.5) kg.",
      source_name = "Weight"
    ),
    BSA = list(
      description = "Body surface area (Du Bois)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained. Table 1 median (IQR): ECMO 1.67 (1.61-1.75), non-ECMO 1.63 (1.51-1.72) m^2.",
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained. Table 1 median (IQR): ECMO 3.00 (2.83-3.20), non-ECMO 2.55 (2.30-2.98) g/dL (source units; 1 g/dL = 10 g/L).",
      source_name = "Albumin"
    ),
    TPRO = list(
      description = "Serum total protein concentration",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as 'serum protein level' (Methods 2.5); not retained. Table 1 median (IQR): ECMO 5.15 (4.88-5.75), non-ECMO 5.05 (4.70-5.75) g/dL (source units; 1 g/dL = 10 g/L).",
      source_name = "Protein"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained directly (it enters only through the CKD-EPI eGFR). Table 1 median (IQR): ECMO 0.820 (0.518-1.15), non-ECMO 0.615 (0.458-1.43) mg/dL.",
      source_name = "Scr"
    ),
    CYSC = list(
      description = "Serum cystatin C concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5); not retained. Table 1 prints the unit as mg/dL (medians 1.48 and 1.34), which is numerically consistent with the usual mg/L scale for cystatin C; the printed label is kept here as published.",
      source_name = "Cystatin C"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26,
    n_studies = 1,
    n_observations = 125,
    age_range = "Adults; median (IQR) 64.0 (56.3-66.5) years on ECMO and 72.0 (66.0-80.3) years without ECMO (Table 1).",
    weight_range = "Median (IQR) 63.5 (61.9-66.3) kg on ECMO and 54.4 (50.5-64.5) kg without ECMO (Table 1).",
    sex_female_pct = "31% (8 of 26: 4 of 8 on ECMO, 4 of 18 without ECMO).",
    race_ethnicity = "Korean (single-centre study in South Korea).",
    disease_state = paste(
      "Critically ill ICU adults receiving meropenem for empirical",
      "management of sepsis of unknown source, nosocomial infection, or",
      "prophylaxis during ECMO. 8 on ECMO (7 veno-arterial, 1 veno-venous);",
      "continuous renal replacement therapy in 1 of 18 without ECMO and in",
      "1 of 8 (Results 3.1) or 3 of 8 (Table 1) on ECMO. APACHE II median",
      "21.0 (ECMO) vs 16.0; SOFA 9.50 vs 5.00."
    ),
    renal_function = "CKD-EPI eGFR median (IQR) 87.7 (70.0-105) mL/min/1.73 m^2 on ECMO and 91.6 (45.6-103) without ECMO; Cockcroft-Gault CLCR 76.9 (59.5-105) and 73.4 (32.7-92.6) mL/min (Table 1).",
    dose_range = "500 or 1000 mg meropenem as a 30-min IV infusion every 8 or 12 h. Simulated regimens: 0.5, 1 and 2 g every 8 or 12 h over 0.5, 1, 2 or 3 h.",
    regions = "South Korea (Hallym University Sacred Heart Hospital, Anyang), September 2020 to April 2021.",
    notes = paste(
      "Five samples after the first dose following enrolment and two at",
      "steady state after the fourth or fifth dose; 125 samples built the",
      "model and 44 trough/peak samples externally validated it. Plasma",
      "meropenem by HPLC-MS/MS (LLOQ 0.2 mg/L). NONMEM 7.5 FOCE-I; PsN",
      "5.2.6 for the stepwise covariate search, prediction- and",
      "variability-corrected VPC and a 2000-sample bootstrap. OFV of one-,",
      "two- and three-compartment base models printed as 689.840, 6540.693",
      "and 640.694 (the two-compartment value is evidently 640.693); final",
      "model OFV 611.402."
    )
  )

  ini({
    # ---- Structural parameters (Lee 2021 Table 2, final model) ----
    # Reported as linear-scale typical values; log-transformed here.
    lcl <- log(6.37); label("Clearance at the reference CKD-EPI eGFR of 91.57 mL/min/1.73 m^2 (L/h)") # Table 2 theta1 = 6.37 L/h (RSE 7.41%; bootstrap 6.32, 95% CI 5.42-7.23)
    lvc <- log(9.07); label("Central volume of distribution (L)") # Table 2 VC = 9.07 L (RSE 12.2%; bootstrap 8.97, 95% CI 3.92-12.0)
    lq <- log(10.7); label("Intercompartmental clearance (L/h)") # Table 2 Q = 10.7 L/h (RSE 21.5%; bootstrap 10.6, 95% CI 4.73-31.0)
    lvp <- log(7.91); label("Peripheral volume of distribution (L)") # Table 2 VP = 7.91 L (RSE 13.6%; bootstrap 8.17, 95% CI 5.35-11.1)

    # ---- Covariate effect (Lee 2021 Table 2 structural-model row) ----
    # CL = theta1 * (1 + theta2 * (CE - 91.57)), CE = CKD-EPI eGFR. CL stays
    # positive for any CE above 91.57 - 1/0.00925 = -16.5, i.e. for every
    # physiological eGFR.
    e_crcl_cl <- 0.00925; label("Linear slope of CKD-EPI eGFR on CL (per mL/min/1.73 m^2)") # Table 2 theta2 = 0.00925 (RSE 10.3%; bootstrap 0.00932, 95% CI 0.00680-0.0110)

    # ---- Plasma protein binding ----
    # Methods 2.7: 'The parameter f was fixed at 98%.' Not estimated from the
    # PK data, hence fixed(). fu converts the total plasma concentration Cc
    # to the free concentration Cu that drives the fT>MIC targets.
    fu <- fixed(0.98); label("Fraction unbound in plasma (unitless)") # Methods 2.7, PD Target Attainment

    # ---- Inter-individual variability (Lee 2021 Table 2) ----
    # Methods 2.5: theta_i = theta * exp(eta_i), eta ~ N(0, omega^2). Table 2
    # reports the IIV as percentages, converted to the log-scale variance via
    # omega^2 = log(CV^2 + 1). No IIV on Q; no covariance reported.
    etalcl ~ 0.0940330 # Table 2: IIV CL = 31.4% (RSE 15.8%, shrinkage 3.70%; bootstrap 29.9, 95% CI 18.0-38.7) -> log(1 + 0.314^2)
    etalvc ~ 0.1740340 # Table 2: IIV VC = 43.6% (RSE 22.5%, shrinkage 14.7%; bootstrap 41.0, 95% CI 0.000-95.4) -> log(1 + 0.436^2)
    etalvp ~ 0.1257124 # Table 2: IIV VP = 36.6% (RSE 21.0%, shrinkage 41.3%; bootstrap 34.5, 95% CI 0.000-55.7) -> log(1 + 0.366^2)

    # ---- Residual error (Lee 2021 Table 2) ----
    # Methods 2.5: 'A power parameter was tested to allow for nonlinear
    # heteroscedastic variances'; Results 3.2 retains it. Encoded as
    # 'Cc ~ pow(propSd, powExp)', i.e. SD(Cobs - Cpred) = propSd * Cpred^powExp.
    propSd <- 0.246; label("Power-error SD coefficient (fraction)") # Table 2: proportional error 24.6% (RSE 29.3%, shrinkage 24.2%; bootstrap 24.1, 95% CI 10.3-41.5)
    powExp <- 0.865; label("Power-error exponent (unitless)") # Table 2: power parameter 0.865 (RSE 10.0%; bootstrap 0.897, 95% CI 0.533-1.38)
  })
  model({
    # 1. Individual PK parameters. CKD-EPI eGFR acts linearly on CL only.
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 91.57))
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    # 2. ODE system, written with flows and volumes (NONMEM ADVAN3 TRANS4
    #    parameterization). Meropenem is given by IV infusion only, so doses
    #    go straight to `central`.
    d/dt(central) <- -cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # 3. Observations. Cc is total plasma meropenem (mg/L; dose in mg, volume
    #    in L). Cu is the free concentration compared against the MIC; it is a
    #    deterministic transform of Cc and carries no residual error.
    Cc <- central / vc
    Cu <- fu * Cc

    Cc ~ pow(propSd, powExp)
  })
}
