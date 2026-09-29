Zhang_2020_teicoplanin <- function() {
  description <- "One-compartment IV population PK model for teicoplanin in 159 hospitalised Chinese children aged 1 month to 14 years (Zhang 2020). Clearance carries a linear body-weight term and an exponential serum-creatinine term, CL = 0.0694 * (1 + 2.82 * WT/16.71) * 0.882^(SCr/29.075), and the volume of distribution an exponential body-weight term, V = 1.39 * 1.75^(WT/16.71); 16.71 kg and 29.075 umol/L are the cohort means. For a child at those means CL = 0.234 L/h (0.014 L/h/kg) and V = 2.43 L (0.15 L/kg). Interindividual variability is exponential on CL (65.9% CV) and V (61.0% CV), and residual variability is 7.0% proportional."
  reference <- "Zhang T, Sun D, Shu Z, Duan Z, Liu Y, Du Q, Zhang Y, Dong Y, Wang T, Hu S, Cheng H, Dong Y. Population Pharmacokinetics and Model-Based Dosing Optimization of Teicoplanin in Pediatric Patients. Front Pharmacol. 2020;11:594562. doi:10.3389/fphar.2020.594562"
  vignette <- "Zhang_2020_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological matrix.
  # Zhang 2020 Methods: teicoplanin in serum by a validated HPLC method
  # (calibration 2.5-100 mg/L, LLOQ 2.5 mg/L).
  compartmentData <- list(
    central = list(analyte = "teicoplanin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Enters CL linearly as (1 + 2.82 * WT/16.71) and V exponentially as 1.75^(WT/16.71) (Zhang 2020 Results, final-model equations). 16.71 kg is the model-building cohort mean (Table 1: 16.7 +/- 10.1 kg, median 14.8, range 2.9-69.0). Both terms are anchored at WT = 0 rather than at the cohort mean, so the typical values lcl and lvc are NOT the values for a typical child; see ini().",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Enters CL exponentially as 0.882^(SCr/29.075) (Zhang 2020 Results, final-model equation). 29.075 umol/L is the model-building cohort mean (Table 1: 29.1 +/- 17.3 umol/L, median 26.0, range 10.0-139.0). Higher creatinine lowers clearance. If no creatinine reading fell within +/- 48 h of dosing the closest available reading was imputed (12 of 236 samples, 5.1%).",
      source_name = "SCr"
    )
  )

  # Screened in the stepwise covariate search (Zhang 2020 Methods, Supplementary
  # Table S3) but NOT retained in the final model. Documentation only -- none of
  # these appears in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Table 1: 4.1 +/- 3.4 years (median 3.7, range 0.2-14.0). Age on V entered the full model in forward selection (Supplementary Table S3 model 3, dOFV -28.6) but was removed in backward elimination (model 7, dOFV +0.014)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Table 1: 72 of 159 (45.3%) female. Screened but not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      notes = "Table 1: 87.8 +/- 47.2 mL/min. Computed with the Cockcroft-Gault formula because height was unavailable for most children; screened but not retained. The Discussion attributes this to Cockcroft-Gault overestimating renal function in small children."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "(not reported)",
      type = "continuous",
      notes = "Screened (Methods) but not retained; no summary statistics reported."
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "(not reported)",
      type = "continuous",
      notes = "Screened (Methods) but not retained; no summary statistics reported."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "(not reported)",
      type = "continuous",
      notes = "Screened (Methods) but not retained; no summary statistics reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 159L,
    n_studies = 1L,
    n_concentrations = 236L,
    age_range = "Mean 4.1 +/- 3.4 years, median 3.7, range 0.2-14.0 years; 51 (32.1%) under 2 years, 98 (61.6%) 2-10 years, 10 (6.3%) 10 years and over (Zhang 2020 Table 1). Eligibility was 1 month to 18 years.",
    weight_range = "Mean 16.7 +/- 10.1 kg, median 14.8, range 2.9-69.0 kg (Table 1)",
    sex_female_pct = 45.3,
    race_ethnicity = "Chinese (two hospitals in Xi'an, China); described by the authors as an Asian paediatric population",
    disease_state = "Hospitalised children receiving teicoplanin for proven or suspected MRSA infection. Indications (not mutually exclusive): respiratory tract infection 97.5%, sepsis 24.5%, bacteraemia 12.6%, bone and joint infection 6.9%. Comorbidities: malignant haematological disease 57.2%, congenital heart disease 15.1%, myocardial injury 13.8%. 30.2% ventilated, 25.2% admitted to intensive care (Table 1).",
    dose_range = "Three loading doses of 10 mg/kg every 12 h, then 6-10 mg/kg/day (label regimen). Received: loading dose 9.8 +/- 1.4 mg/kg (range 5.2-16.0); daily maintenance dose 9.5 +/- 1.2 mg/kg (range 5.2-12.9) (Table 1). Infusion duration not reported.",
    regions = "China (First Affiliated Hospital and Affiliated Children Hospital of Xi'an Jiaotong University)",
    renal_function = "Serum creatinine 29.1 +/- 17.3 umol/L (median 26.0, range 10.0-139.0); Cockcroft-Gault creatinine clearance 87.8 +/- 47.2 mL/min (Table 1).",
    co_medication = "Other antibacterials: ceftriaxone 42.8%, meropenem 34.0%, imipenem-cilastatin 45.3%, cefoperazone-sulbactam 19.5%; loop diuretic 42.8% (Table 1). Nephrotoxic co-medication was screened but not retained.",
    notes = "Retrospective study, March 2017 to November 2019. Sparse data: 236 concentrations, 1.5 per child, 212 (89.8%) of them TDM troughs drawn within 30 min before a dose at steady state; six below the 2.5 mg/L LLOQ were set to 2.5 mg/L. NONMEM 7.2, ADVAN1 TRANS2, FOCE-I. Final-model OFV 971.014. Shrinkage 26.9% (CL), 19.8% (V), 24.4% (residual). The model was externally evaluated on a separate cohort of 66 children (89 concentrations) from the same hospitals; that cohort is not part of the fit."
  )

  ini({
    # Structural parameters from Zhang 2020 Table 2 ('Final model' column),
    # confirmed by the final-model equations printed in the Results:
    #   CL (L/h) = 0.0694 x (1 + theta1 x WT/16.71) x theta2^(SCr/29.075) x e^eta1
    #   Vd (L)   = 1.39 x theta3^(WT/16.71) x e^eta2
    # Both are anchored at WT = 0 (and SCr = 0), so these two values are not the
    # values of any real child. The abstract's 'clearance ... 0.694 L/h' drops a
    # zero: Table 2 and the equation both give 0.0694, and only 0.0694
    # reproduces the abstract's own '0.784 L/h/70 kg' (see e_creat_cl below).
    lcl <- log(0.0694); label("Clearance at WT = 0 and SCr = 0 (L/h)")                 # Table 2: CL = 0.0694 L/h (RSE 11.3%; bootstrap mean 0.0718, 95% CI 0.0453-0.0983)
    lvc <- log(1.39); label("Volume of distribution at WT = 0 (L)")                  # Table 2: Vd = 1.39 L (RSE 11.0%; bootstrap mean 1.77, 95% CI 1.34-2.20)

    # Covariate effects. The printed equation image shows theta2 and theta3 as
    # BASES raised to the covariate ratio (theta^(cov/ref)), not as power
    # exponents on the ratio ((cov/ref)^theta). Three of the paper's own derived
    # numbers reproduce only under the base reading:
    #   CL at WT 16.71, SCr 29.075 = 0.0694 * 3.82 * 0.882 = 0.2338 L/h
    #     = 0.0140 L/h/kg (Discussion: '0.014 L/h/kg'); the power reading gives 0.0159.
    #   CL at WT 70, SCr 29.075 = 0.0694 * 12.813 * 0.882 = 0.784 L/h
    #     (Abstract: '0.784 L/h/70 kg'); the power reading gives 0.889.
    #   Vd at WT 16.71 = 1.39 * 1.75 = 2.43 L = 0.146 L/kg
    #     (Discussion: 'Vd in this study (0.15 L/kg)'); the power reading gives 0.083.
    e_wt_cl <- 2.82; label("Linear slope of CL on WT/16.71 (unitless)")               # Table 2: theta_wt on CL = 2.82 (RSE 20.6%; bootstrap mean 3.62, 95% CI 1.21-6.03)
    e_creat_cl <- 0.882; label("Base of the exponential SCr effect on CL, theta^(SCr/29.075) (unitless)") # Table 2: theta_SCr on CL = 0.882 (RSE 5.0%; bootstrap mean 0.794, 95% CI 0.688-0.9)
    e_wt_vc <- 1.75; label("Base of the exponential WT effect on Vd, theta^(WT/16.71) (unitless)")      # Table 2: theta_wt on Vd = 1.75 (RSE 6.3%; bootstrap mean 1.76, 95% CI 1.29-2.23)

    # IIV: exponential model (Methods), reported in Table 2 as CV%. Converted with
    # omega^2 = log(1 + CV^2); the paper gives no CI, equation exponent or
    # variance column that would settle the convention otherwise.
    etalcl ~ 0.3607 # Table 2: CV-CL = 65.9% (RSE 17.6%); log(1 + 0.659^2) = 0.3607
    etalvc ~ 0.3163 # Table 2: CV-Vd = 61.0% (RSE 42.5%); log(1 + 0.610^2) = 0.3163

    # Residual error: Table 2 reports 'Residual variability (%)', 'CV-sigma' =
    # 7.0%. The Results say an additive model was selected, and a CV% label on an
    # additive error means additive on the log-transformed scale, i.e.
    # proportional on the linear scale. A 7.0 mg/L additive SD is ruled out by
    # Figure 3: its lower 5th-percentile band stays at 0-5 mg/L, whereas a
    # 7 mg/L additive SD at a median of ~10 mg/L would put it well below zero.
    propSd <- 0.07; label("Proportional residual error (fraction)")                    # Table 2: CV-sigma = 7.0% (RSE 21.9%; bootstrap mean 8.5, 95% CI 5.1-11.9)
  })
  model({
    # Individual PK parameters (Zhang 2020 Results, final-model equations).
    # 16.71 kg and 29.075 umol/L are the model-building cohort means (Table 1).
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * WT / 16.71) * e_creat_cl^(CREAT / 29.075)
    vc <- exp(lvc + etalvc) * e_wt_vc^(WT / 16.71)

    kel <- cl / vc

    # One-compartment model with first-order elimination (ADVAN1 TRANS2).
    # Doses (mg) enter central as IV infusions; the infusion duration was not
    # reported and is set in the event table.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
