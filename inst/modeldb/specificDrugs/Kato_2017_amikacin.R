Kato_2017_amikacin <- function() {
  description <- "Two-compartment IV population PK model for amikacin in hospitalized Japanese adults treated for Gram-negative and other infections, with clearance directly proportional to Cockcroft-Gault creatinine clearance and central volume directly proportional to total body weight (Kato 2017, Phoenix NLME, n = 35 patients, 65 therapeutic drug monitoring samples)."
  reference <- "Kato H, Hagihara M, Hirai J, Sakanashi D, Suematsu H, Nishiyama N, Koizumi Y, Yamagishi Y, Matsuura K, Mikamo H. Evaluation of Amikacin Pharmacokinetics and Pharmacodynamics for Optimal Initial Dosing Regimen. Drugs R D. 2017;17(1):177-187. doi:10.1007/s40268-016-0165-5"
  vignette <- "Kato_2017_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Source name CCr. Estimated by the Cockcroft-Gault equation (Kato 2017 Data Collection; Table 1 footnote a). Taken at the initiation of treatment (Covariate Model Building), so time-fixed per subject. Enters clearance as the through-origin linear ratio CL = theta1 * (CCr / 60.0) (Table 3 footnote a; Table 4 footnote), where 60.0 mL/min is the Table 1 cohort median; clearance tends to zero as CRCL tends to zero, so CRCL must be strictly positive. Observed range 6.1-144 mL/min.",
      source_name = "CCr"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source name 'weight' (total body weight; ideal and lean body weight were screened separately and not retained). Enters the central volume as the through-origin linear ratio V1 = theta2 * (weight / 49.5) (Table 3 footnote a; Table 4 footnote), where 49.5 kg is the Table 1 cohort median. Observed range 32.5-78 kg.",
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened (gender) in the stepwise covariate search but not retained (Kato 2017 Covariate Model Building)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    n_samples = 65L,
    age_range = "15-95 years",
    age_median = "70.0 years (mean 67.5 +/- 17.3)",
    weight_range = "32.5-78 kg",
    weight_median = "49.5 kg (mean 50.4 +/- 10.4)",
    sex_female_pct = 28.6,
    race_ethnicity = "Japanese (single Japanese university hospital); not tabulated",
    disease_state = "Inpatients older than 15 years treated with amikacin for at least 3 days: urinary tract infection 8, pneumonia 10, bacteremia 11, pneumonia plus bacteremia 2, peritoneum inflammation 2, febrile neutropenia 2. Pseudomonas aeruginosa was the commonest isolate (48.5% of 33 isolates). Almost all patients (97.1%) received combination antibiotic therapy. Baseline albumin 2.7 +/- 0.6 g/dL.",
    renal_function = "Cockcroft-Gault creatinine clearance mean 60.3 +/- 33.9 mL/min, median 60.0 (range 6.1-144); patients on continuous renal replacement therapy at the onset of amikacin therapy were excluded.",
    dose_range = "Amikacin 200-1000 mg/day IV once daily (median 400 mg/day; 9.6 mg/kg/day, range 3.1-17.7), infused over 0.5 h (30 patients) or 1 h (5 patients).",
    regions = "Japan (Aichi Medical University Hospital)",
    notes = "Single-centre retrospective therapeutic drug monitoring study, September 2009 - December 2014 (Kato 2017 Patients; Table 1). One or two serum samples per patient: a trough immediately before a dose and a peak 1 h after the start of infusion, mostly within 3 days of the start of therapy. Amikacin assayed by fluorescence polarization immunoassay (limit of quantification 0.8 mg/L, calibration range 0.8-40 mg/L, intra- and inter-assay CV < 6%)."
  )

  ini({
    # Structural parameters: Kato 2017 Table 4 'Final estimates' and its
    # footnote 'Final model: CL (L/h) = theta1(CCr/60.0), V1 (L) =
    # theta2(weight/49.5), V2 (L) = 7.72, Q (L/h) = 7.0'.
    lcl <- log(2.41)
    label("Clearance theta1 at CRCL = 60 mL/min (L/h)")
    # Kato 2017 Table 4: theta1 = 2.41 L/h (bootstrap median 2.50, 95% CI 1.97-3.31)
    lvc <- log(10.70)
    label("Central volume V1 theta2 at WT = 49.5 kg (L)")
    # Kato 2017 Table 4: theta2 = 10.70 L (bootstrap median 10.26, 95% CI 5.99-14.69)
    lvp <- log(7.72)
    label("Peripheral volume V2 (L)")
    # Kato 2017 Table 4 footnote: V2 (L) = 7.72 (no uncertainty or bootstrap reported)
    lq <- log(7.0)
    label("Intercompartmental clearance Q (L/h)")
    # Kato 2017 Table 4 footnote: Q (L/h) = 7.0 (no uncertainty or bootstrap reported)

    # Between-subject variability: exponential (log-normal) random effects on
    # CL and V1 only; IIV on Q and V2 'did not improve model fitting' and was
    # not included (Basic Model). Table 4 reports omega as a percentage,
    # converted here by omega^2 = log(CV^2 + 1). No correlation is reported.
    etalcl ~ 0.229577 # Kato 2017 Table 4 omega(CL) 50.8%: log(0.508^2 + 1)
    etalvc ~ 0.080750 # Kato 2017 Table 4 omega(V1) 29.0%: log(0.290^2 + 1)

    # Residual error: multiplicative, 'Cp = F*(1 + epsilon)' (Basic Model),
    # with sigma = 0.25 from Table 4. Table 4 prints the unit as ug/mL, but the
    # printed error equation is proportional and Figure 1D (IWRES within +/-2
    # at IPRED up to ~65 mg/L) is only consistent with a 25% proportional SD;
    # see the vignette 'Assumptions and deviations'.
    propSd <- 0.25
    label("Proportional residual SD (fraction)")
    # Kato 2017 Table 4: sigma = 0.25 (bootstrap median 0.23, 95% CI 0.01-0.33)
  })

  model({
    # Reference covariate values (Table 3 footnote a): the Table 1 cohort
    # medians of creatinine clearance (60.0 mL/min) and body weight (49.5 kg).
    crcl_ref <- 60.0
    wt_ref <- 49.5

    # Individual parameters. CL and V1 are directly proportional to their
    # covariates (no exponent and no intercept in the printed equations).
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order elimination from the central
    # compartment; amikacin is given by IV infusion into the central
    # compartment (Population PK Parameters of Amikacin).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> Cc in mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
