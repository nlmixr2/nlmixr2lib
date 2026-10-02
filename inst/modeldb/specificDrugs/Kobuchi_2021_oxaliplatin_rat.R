Kobuchi_2021_oxaliplatin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model with linear",
    "elimination for intact (unbiotransformed) oxaliplatin in plasma after",
    "a single intravenous bolus of 3 or 8 mg/kg to male Wistar rats with",
    "normal renal function or mild / severe acute kidney injury induced by",
    "30 / 60 min renal ischemia-reperfusion (Kobuchi 2021). Parameters are",
    "per kg body weight (dose in mg/kg), exponential IIV on the central",
    "volume, clearance and intercompartmental clearance, proportional",
    "residual error. The final model carries no covariates: the paper's",
    "renal-function simulation (Figure 5) substituted a post hoc",
    "clearance-versus-plasma-creatinine regression whose coefficients are",
    "not printed (see the vignette)."
  )
  reference <- paste(
    "Kobuchi S, Kai M, Ito Y. Population Pharmacokinetic Model-Based",
    "Evaluation of Intact Oxaliplatin in Rats with Acute Kidney Injury.",
    "Cancers (Basel). 2021;13(24):6382. doi:10.3390/cancers13246382",
    "(PMCID PMC8699120).",
    sep = " "
  )
  vignette <- "Kobuchi_2021_oxaliplatin_rat"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Plasma creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not part of the estimated model (Table 3 lists no covariate",
        "effect). Section 2.5: the Figure 5 simulations set CL from 'the",
        "regression equations of Cr level and post hoc CL' over plasma Cr",
        "0.3-2.5 mg/dL, but the regression coefficients are not printed.",
        "Group means (Table 1): 0.27 (normal), 0.54 (mild), 0.95 (severe)",
        "mg/dL. Blood urea nitrogen and creatinine clearance were correlated",
        "with NCA CLtot (Figure 3) but were likewise not used in the model."
      ),
      source_name = "Cr"
    )
  )

  compartmentData <- list(
    central = list(analyte = "oxaliplatin (intact)", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "oxaliplatin (intact)", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 30,
    n_studies = 1,
    sex = "male",
    age_range = "10 weeks at purchase",
    weight_range = "approximately 300 g (Section 3.2 illustrates the 8 mg/kg dose as 2400 ug for a 300 g rat)",
    sex_female_pct = 0,
    disease_state = paste(
      "Normal renal function, or mild (30 min) / severe (60 min) bilateral renal",
      "ischemia followed by 24 h reperfusion (acute kidney injury)"
    ),
    dose_range = "Oxaliplatin (Elplat) 3 or 8 mg/kg single intravenous bolus into the jugular vein",
    regions = "Japan",
    renal_function = paste(
      "Table 1 plasma creatinine (mean +/- SD): normal 0.27 +/- 0.02, mild",
      "0.54 +/- 0.21, severe 0.95 +/- 0.17 mg/dL; creatinine clearance 4.2 +/-",
      "1.3, 2.2 +/- 0.7, 1.6 +/- 0.7 mL/min/kg."
    ),
    notes = paste(
      "Section 2.3: 30 rats split into normal, mild and severe renal-failure",
      "groups, each dosed at 3 or 8 mg/kg (n = 5 per dose group). Plasma sampled",
      "at 3, 5, 10, 20, 30, 45 min and 1, 1.5, 2 h; intact oxaliplatin measured",
      "by LC-MS/MS. Model fitted in Phoenix NLME 8.2 (FOCE-ELS)."
    )
  )

  ini({
    # Table 3 'Fixed effect parameters (theta)'. Parameters are per kg body
    # weight; the dose is supplied in mg/kg.
    lvc <- log(0.44); label("Central volume V (L/kg)") # Table 3: V = 0.44 L/kg (CV% 7.4; bootstrap median 0.44, 0.39-0.49)
    lvp <- log(2.26); label("Peripheral volume V2 (L/kg)") # Table 3: V2 = 2.26 L/kg (CV% 19.1; bootstrap median 2.28, 1.62-3.16)
    lcl <- log(1.76); label("Clearance from the central compartment CL (L/h/kg)") # Table 3: CL = 1.76 L/h/kg (CV% 8.8; bootstrap median 1.76, 1.50-2.03)
    lq <- log(1.0); label("Intercompartmental clearance CL2 (L/h/kg)") # Table 3: CL2 = 1.0 L/h/kg (CV% 10.0; bootstrap median 1.0, 0.84-1.10)

    # Table 3 'Inter-individual variability (omega)', reported in %; the
    # exponential IIV model is stated in Section 2.4. The percentages are read
    # as omega x 100 (omega^2 = (pct/100)^2), the convention the same group's
    # Phoenix NLME tables were shown to use in Kobuchi 2025 dapagliflozin; the
    # lognormal-CV reading, omega^2 = log(1 + (pct/100)^2), differs by under 7%
    # in variance (see vignette).
    etalvc ~ 0.140625 # Table 3: IIV V = 37.5% (CV% 19.3; bootstrap 29.6-43.9); 0.375^2
    etalcl ~ 0.093025 # Table 3: IIV CL = 30.5% (CV% 26.0; bootstrap 21.6-37.5); 0.305^2
    etalq ~ 0.099225 # Table 3: IIV CL2 = 31.5% (CV% 18.7; bootstrap 24.8-36.6); 0.315^2

    propSd <- 0.149; label("Proportional residual error (fraction)") # Table 3 'Residual variability (sigma)': C = 14.9% (CV% 6.9; bootstrap 13.1-17.0); proportional per Section 2.4
  })

  model({
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
