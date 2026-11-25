Kawata_2019_gefitinib <- function() {
  description <- "One-compartment population PK model for oral gefitinib in Japanese adults with advanced or recurrent non-small-cell lung cancer (Kawata 2019). First-order absorption with lag time; power-form effects of serum alpha1-acid glycoprotein (reference 121 mg/dL) and body weight (reference 55.5 kg) on CL/F and V/F, multiplicative factors on CL/F for age >= 65 years and concomitant CYP3A4 inducers, and a time-varying multiplicative factor on CL/F at the time of interstitial lung disease (ILD) diagnosis. Exponential IIV on CL/F, V/F, ka and tlag; additive residual error on log-transformed concentrations."
  reference <- "Kawata T, Higashimori M, Itoh Y, Tomkinson H, Johnson MG, Tang W, Nyberg F, Jiang H, Tanigawara Y. Gefitinib exposure and occurrence of interstitial lung disease in Japanese patients with non-small-cell lung cancer. Cancer Chemother Pharmacol. 2019;83(5):849-858. doi:10.1007/s00280-019-03788-4"
  vignette <- "Kawata_2019_gefitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "gefitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gefitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AAG = list(
      description = "Serum alpha1-acid glycoprotein (AGP) concentration.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Kawata 2019 reports AGP in mg/dL (Table 1: 121.3 +/- 52.4 mg/dL, range 31-327) and centres the power-form effects on 121 mg/dL (Table 2 CL/F and V/F equations). The canonical column is in g/L, so the model uses the equivalent reference 1.21 g/L: AAG [g/L] = AGP [mg/dL] / 100. Measured at first dosing and at steady state; not measured at the time of ILD development (Discussion).",
      source_name = "AGP"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effects on CL/F and V/F centred on the cohort mean 55.5 kg (Table 1; Table 2 CL/F and V/F equations).",
      source_name = "WGT"
    ),
    AGE_GE65 = list(
      description = "Age >= 65 years indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (< 65 years)",
      notes = "Table 2 footnote: 'AGE >= 65 is 0 for patients aged < 65 years and 1 for patients aged >= 65 years'. 180 of 336 patients were >= 65 years (Table 1). Multiplicative factor theta3^AGE_GE65 on CL/F.",
      source_name = "AGE >= 65"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant use of a CYP3A4 inducer.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CYP3A4 inducer)",
      notes = "Table 2 footnote: 'IND is 0 for no concomitant use with CYP3A4 inducers and 1 for concomitant use with CYP3A4 inducers'. 23 of 336 patients (Table 1). The paper does not list which inducers or strengths were pooled into IND = 1, and detailed dosing histories of the concomitant drugs were not collected (Discussion). Multiplicative factor theta4^CONMED_CYP3A4_IND on CL/F.",
      source_name = "IND"
    ),
    DIS_ILD_ACUTE = list(
      description = "Record taken at the time of interstitial lung disease (ILD) diagnosis.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any time other than ILD development, including all records of patients who never developed ILD)",
      notes = "Table 2 footnote: 'PERI is 0 at the time of no ILD diagnosis and 1 at the time of ILD diagnosis'. Time-varying within subject: 37 concentrations from 27 of the 51 ILD patients were drawn at ILD development (Methods, Plasma concentrations). A subject-level 'patients with ILD' indicator was also tested and not retained (Online Resource, Population pharmacokinetic model development). Multiplicative factor theta6^DIS_ILD_ACUTE on CL/F.",
      source_name = "PERI"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested as a covariate and not retained in the final model (Online Resource, Population pharmacokinetic model development).",
      source_name = "sex"
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant use of a CYP3A4 inhibitor.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CYP3A4 inhibitor)",
      notes = "Tested and not retained (Online Resource; Discussion notes the inconsistency with the itraconazole interaction study). 24 of 336 patients (Table 1).",
      source_name = "concomitant CYP3A4 inhibitors"
    ),
    CONMED_PPI = list(
      description = "Concomitant use of a proton-pump inhibitor and/or H2-receptor antagonist.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no acid-reducing agent)",
      notes = "Tested as a pooled PPI and/or H2-antagonist indicator and not retained (Online Resource; Discussion). 115 of 336 patients (Table 1).",
      source_name = "concomitant PPIs and/or H2 antagonists"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a power model centred on the mean and not retained (Online Resource). Total protein, AST, ALT, ALP, total bilirubin and creatinine clearance were also screened and not retained.",
      source_name = "serum albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 336L,
    n_studies = 1L,
    n_observations = 1891L,
    age_range = "37-87 years (mean 64.8 +/- 9.1)",
    weight_range = "33-100 kg (mean 55.5 +/- 11.0)",
    sex_female_pct = 45.2,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Advanced or recurrent non-small-cell lung cancer after at least one prior chemotherapy regimen, irrespective of EGFR mutation status; 51 patients developed interstitial lung disease (ILD) and 285 did not.",
    dose_range = "Gefitinib 250 mg orally once daily.",
    regions = "Japan (50 centres, November 2003 to February 2006).",
    co_medication = "CYP3A4 inducers 23 patients, CYP3A4 inhibitors 24, proton-pump inhibitors and/or H2 antagonists 115 (Table 1).",
    notes = "Population PK analysis set of a prospective cohort with a nested ILD case-control study (NCT00252759). Samples at 1-3, 3-8 and 24 h after the first dose and at steady state (days 10-15), plus samples at ILD development where possible: 943 first-dose concentrations from 332 patients, 911 steady-state concentrations from 307 patients, and 37 at ILD development from 27 patients. Serum AGP 121.3 +/- 52.4 mg/dL (31-327). Fitted in NONMEM V (FO method) to log-transformed concentrations (Kawata 2019 Table 1, Methods and Online Resource)."
  )

  ini({
    # Structural parameters -- Kawata 2019 Table 2 final estimates. CL/F and V/F
    # are apparent oral values (bioavailability not identifiable); f(depot) = 1.
    lcl <- log(28.6); label("Apparent oral clearance CL/F for the reference patient (L/h)") # Table 2: theta1 = 28.6 L/h
    lvc <- log(1540); label("Apparent volume of distribution V/F for the reference patient (L)") # Table 2: theta7 = 1540 L
    lka <- log(0.819); label("First-order absorption rate constant ka (1/h)") # Table 2: ka = theta10 = 0.819 1/h
    ltlag <- log(0.863); label("Absorption lag time tlag (h)") # Table 2: tlag = theta11 = 0.863 h

    # Covariate effects on CL/F -- Table 2 equation
    # CL/F = theta1 * (AGP/121)^theta2 * theta3^(AGE>=65) * theta4^IND * (WGT/55.5)^theta5 * theta6^PERI
    e_aag_cl <- -0.502; label("Power exponent of AGP/121 mg/dL on CL/F (unitless)") # Table 2: theta2 = -0.502
    e_age_ge65_cl <- 0.759; label("Multiplicative factor on CL/F for age >= 65 years (unitless)") # Table 2: theta3 = 0.759
    e_cyp3a4_ind_cl <- 1.82; label("Multiplicative factor on CL/F for concomitant CYP3A4 inducer (unitless)") # Table 2: theta4 = 1.82
    e_wt_cl <- 0.503; label("Power exponent of WT/55.5 kg on CL/F (unitless)") # Table 2: theta5 = 0.503
    e_dis_ild_acute_cl <- 0.743; label("Multiplicative factor on CL/F at the time of ILD diagnosis (unitless)") # Table 2: theta6 = 0.743

    # Covariate effects on V/F -- Table 2 equation
    # V/F = theta7 * (AGP/121)^theta8 * (WGT/55.5)^theta9
    e_aag_vc <- -0.447; label("Power exponent of AGP/121 mg/dL on V/F (unitless)") # Table 2: theta8 = -0.447
    e_wt_vc <- 0.766; label("Power exponent of WT/55.5 kg on V/F (unitless)") # Table 2: theta9 = 0.766

    # IIV -- exponential model, no covariances (Online Resource). Table 2 prints
    # each row as a percentage, read as a log-normal CV:
    # variance = log(1 + CV^2).
    etalcl ~ 0.217567 # Table 2: omega CL/F = 49.3% -> log(1 + 0.493^2)
    etalvc ~ 0.138217 # Table 2: omega V/F = 38.5% -> log(1 + 0.385^2)
    etalka ~ 2.08885 # Table 2: omega ka = 266% -> log(1 + 2.66^2)
    etaltlag ~ 0.0678689 # Table 2: omega tlag = 26.5% -> log(1 + 0.265^2)

    # Residual error -- additive on log-transformed concentrations (Online
    # Resource), i.e. log-normal on the linear scale. Table 2 sigma = 32.7%,
    # read on the same CV scale as the omegas: SD = sqrt(log(1 + 0.327^2)).
    expSd <- 0.318731; label("Residual SD on the log-concentration scale (unitless)") # Table 2: sigma = 32.7% -> sqrt(log(1 + 0.327^2))
  })

  model({
    # Reference covariate values from the Table 2 equations; AGP 121 mg/dL is
    # 1.21 g/L in the canonical AAG units.
    ref_aag <- 1.21
    ref_wt <- 55.5

    cl <- exp(lcl + etalcl) *
      (AAG / ref_aag)^e_aag_cl *
      e_age_ge65_cl^AGE_GE65 *
      e_cyp3a4_ind_cl^CONMED_CYP3A4_IND *
      (WT / ref_wt)^e_wt_cl *
      e_dis_ild_acute_cl^DIS_ILD_ACUTE
    vc <- exp(lvc + etalvc) * (AAG / ref_aag)^e_aag_vc * (WT / ref_wt)^e_wt_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    # One-compartment model with first-order absorption and lag time
    # (NONMEM ADVAN2 TRANS2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
