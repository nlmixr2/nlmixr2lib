Cojutti_2026_fosfomycin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "fosfomycin trometamol in plasma, joined to prostate tissue by a direct-",
    "response power function (Cprostate = PR * Cc^POW), in 104 men undergoing",
    "endoscopic surgery (HoLEP or TURP) for benign prostatic hyperplasia after a",
    "prophylactic scheme of two 3 g oral doses 12 h apart. Prostate",
    "concentration is algebraic in the plasma concentration; there is no",
    "prostate ODE state. Clearance and volume are apparent (CL/F, V/F); ka and",
    "its IIV are fixed. V/F carries an exponential effect of uncentred body",
    "weight. Constant (additive) residual error on both outputs (Cojutti 2026).",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Berrino PM, Rotaru V, Bianchi L, Schiavina R, Brunocilla E,",
    "Viale P, Pea F. Population Pharmacokinetics of a Novel Oral Fosfomycin",
    "Prophylactic Scheme in the Plasma and Prostate of Patients Undergoing",
    "Endoscopic Surgery for Benign Prostatic Hyperplasia.",
    "Clin Pharmacokinet. 2026;65:109-118.",
    "doi:10.1007/s40262-025-01586-1.",
    sep = " "
  )
  vignette <- "Cojutti_2026_fosfomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds. Prostate tissue is an algebraic observable of
  # the plasma concentration (ESM mlxtran: Cp = PR*(Cc^POW)), not a state;
  # Cprostate is reported in ug/g of tissue.
  compartmentData <- list(
    depot = list(analyte = "fosfomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Applied as an exponential on ABSOLUTE (uncentred) body weight in kg:",
        "V/F = 60.86 * exp(0.02 * BW) (Cojutti 2026 Results, 'Model equations').",
        "The typical value 60.86 L is therefore V/F at 0 kg; at the cohort",
        "median of 78 kg it is 60.86 * exp(0.02 * 78) = 289.6 L, matching the",
        "paper's median individual apparent V of 288.9 L. Methods calls the",
        "covariate model 'a power function', but the printed equation and the",
        "288.9 L median both require the exponential form. Cohort median 78,",
        "range 60-107 kg. Body weight was the only covariate retained (ESM",
        "Table S2: coefficient 0.31, p = 1.42e-3 on V; dOFV 5.49, 32.3%",
        "reduction of the IIV of V).",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and V in the univariate analysis and not retained (ESM Table S2: p = 0.745 on CL, 0.946 on V). Cohort median 70, range 49-84 years."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.339 on CL, 0.353 on V). Cohort median 26.3, range 20.6-41.5 kg/m^2."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.0619 on CL, 0.807 on V)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by CKD-EPI",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.132 on CL, 0.878 on V). Cohort median 78.5, range 34-102 mL/min/1.73 m^2; 102 of 104 patients had eGFR above 50 mL/min, so the data could not inform a renal effect (Discussion)."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.892 on CL, 0.283 on V)."
    ),
    PSA_BL = list(
      description = "Serum prostate-specific antigen before intervention",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.819 on CL, 0.426 on V). Cohort median 3.10, range 0.02-66 ng/mL."
    ),
    PROSTATE_VOL_BL = list(
      description = "Prostate volume at transrectal ultrasound before intervention",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (ESM Table S2: p = 0.577 on CL, 0.746 on V). Cohort median 75, range 20-220 mL."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 104L,
    n_studies = 1L,
    n_observations = 208L,
    age_range = "49-84 years (median 70)",
    weight_range = "60-107 kg (median 78)",
    sex_female_pct = 0,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Men with benign prostatic hyperplasia undergoing endoscopic surgery:",
      "Holmium laser enucleation of the prostate (HoLEP, 85/104, 81.7%) or",
      "transurethral resection of the prostate (TURP, 19/104, 18.3%).",
      "Prostate volume median 75 (20-220) mL; PSA median 3.10 (0.02-66) ng/mL.",
      sep = " "
    ),
    renal_function = "CKD-EPI eGFR median 78.5 (34-102) mL/min/1.73 m^2; 98% (102/104) above 50 mL/min.",
    dose_range = paste(
      "Oral fosfomycin trometamol 3 g twice, 12 h apart, the second dose given",
      "before surgery (3 g at t = -12 h and 3 g at t = 0 h relative to the",
      "end of the prophylactic scheme).",
      sep = " "
    ),
    regions = "Italy (single centre: IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes = paste(
      "One plasma sample and one concomitant prostate biopsy per patient, at a",
      "median of 3.86 (0.42-10.7) h after the second dose; no sampling during",
      "absorption, which is why ka was fixed. Observed median plasma 11.1",
      "(2.1-26.6) mg/L and prostate 11.6 (1.2-603) ug/g; observed prostate-to-",
      "plasma ratio median 1.04 (0.23-5.85). LLOQ 1 mg/L in plasma and",
      "0.1 ug/g in prostate (LC-MS/MS). Fitted in Monolix 2023R1 (SAEM).",
      sep = " "
    )
  )

  ini({
    # Plasma PK -- Table 1 'Fixed effects' (identical to ESM Table S3, 'Model
    # with Weight on V'). CL and V are apparent (oral data only; F = 1).
    lka <- fixed(log(0.92)); label("Absorption rate constant ka (1/h)") # Table 1 ka = 0.92 1/h, fixed (Methods 2.2: value from healthy volunteers, ref 23)
    lcl <- log(14.68); label("Apparent clearance CL/F (L/h)") # Table 1 CL = 14.68 L/h (RSE 14.4%)
    lvc <- log(60.86); label("Apparent volume V/F at 0 kg body weight (L); see e_wt_vc") # Table 1 V = 60.86 L (RSE 44.1%); Results equation V = 60.86 x e^(0.02 x BW)
    e_wt_vc <- 0.02; label("Exponential effect of absolute body weight on V/F (1/kg)") # Table 1 beta_V_weight = 0.02 (RSE 26.9%)

    # Plasma-to-prostate power link -- Table 1; Methods 2.2 model (2)
    # Cpr(t) = PR x Cpl(t)^POW; Results equation Cp = 0.26 x Cc^1.57.
    lrprostate <- log(0.26); label("Prostate / plasma penetration ratio at a plasma concentration of 1 mg/L ((ug/g)/(mg/L)^pow_prostate)") # Table 1 PR = 0.26 (RSE 7.09%); footnote 'prostate/blood penetration ratio at a plasma concentration of 1 mg/L'
    lpow_prostate <- log(1.57); label("Power exponent of plasma concentration in the plasma-to-prostate link (unitless)") # Table 1 POW = 1.57 (RSE 1.80%)

    # Between-subject variability -- Table 1 'Standard deviation of the random
    # effects'. Monolix reports SDs on the log scale; the ESM Table S3 CV column
    # confirms it (omega_V 0.38 -> 39.52% = sqrt(exp(0.38^2) - 1)). nlmixr2
    # takes variances, so each SD is squared. No correlations were retained.
    etalka ~ fixed(1.09^2) # Table 1 omega_ka = 1.09 (SD, no RSE); ESM Table S3 CV 151.02%
    etalcl ~ 0.76^2 # Table 1 omega_CL = 0.76 (SD, RSE 12.4%); ESM Table S3 CV 87.86%
    etalvc ~ 0.38^2 # Table 1 omega_V = 0.38 (SD, RSE 17.8%); ESM Table S3 CV 39.52%
    etalrprostate ~ 0.42^2 # Table 1 omega_PR = 0.42 (SD, RSE 11.8%); ESM Table S3 CV 44.27%
    etalpow_prostate ~ 0.095^2 # Table 1 omega_POW = 0.095 (SD, RSE 15.4%); ESM Table S3 CV 9.5%

    # Residual error -- Table 1 'Error model parameters'; footnote: 'The letter
    # a for error model identifies the estimated parameter for the constant
    # error model'.
    addSd <- 1.35; label("Plasma additive residual error (mg/L)") # Table 1 a_plasma = 1.35 mg/L (RSE 20.5%)
    addSd_Cprostate <- 1.01; label("Prostate additive residual error (ug/g)") # Table 1 a_prostate = 1.01 ug/g (RSE 34.0%)
  })

  model({
    # Individual parameters (exponential random effects, Methods 2.2)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * exp(e_wt_vc * WT)
    rprostate <- exp(lrprostate + etalrprostate)
    pow_prostate <- exp(lpow_prostate + etalpow_prostate)

    kel <- cl / vc

    # ESM mlxtran: Cc = pkmodel(ka, V, Cl), one compartment, first-order
    # absorption and elimination
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Plasma concentration (mg/L) and the ESM mlxtran prostate link
    # Cp = PR*(Cc^POW) (ug/g)
    Cc <- central / vc
    Cprostate <- rprostate * Cc^pow_prostate

    Cc ~ add(addSd)
    Cprostate ~ add(addSd_Cprostate)
  })
}
