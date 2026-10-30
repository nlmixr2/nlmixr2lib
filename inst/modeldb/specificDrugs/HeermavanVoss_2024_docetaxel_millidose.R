HeermavanVoss_2024_docetaxel_millidose <- function() {
  description <- paste(
    "Two-compartment population PK model for a sub-pharmacological",
    "intravenous docetaxel test dose (millidose, 1000 ug bolus) in 29 men",
    "with prostate cancer, fitted separately from the therapeutic-dose data",
    "of the same patients. Clearance and inter-compartmental clearance scale",
    "with (WT/70)^0.75 and both volumes with (WT/70)^1 (fixed exponents).",
    "Inter-individual variability on CL and the peripheral volume;",
    "proportional residual error."
  )
  reference <- paste(
    "Heerma van Voss MR, Notohardjo J, van Dodewaard-de Jong J,",
    "Bloemendal HJ, ter Heine R. A sub-pharmacological test dose does not",
    "predict individual docetaxel exposure in prostate cancer patients.",
    "Cancer Chemother Pharmacol. 2024;94(3):437.",
    "doi:10.1007/s00280-024-04684-2."
  )
  vignette <- "HeermavanVoss_2024_docetaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL and Q with exponent 0.75 and of Vc and Vp",
        "with exponent 1, all referenced to 70 kg (Methods: 'Flow and volume",
        "parameters were allometrically scaled to a total body weight of",
        "70 kg, with allometric coefficients of 0.75 and 1'). Table 1",
        "cohort weight mean 85.5 kg (95% CI 69.4-100.8)."
      ),
      source_name = "total body weight"
    )
  )

  covariatesDataExcluded <- list(
    AAG = list(
      description = "Serum alpha-1-acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested as a covariate on clearance but not retained: 'AAG",
        "concentrations were not significantly correlated with clearance of",
        "docetaxel, and this covariate was, therefore, not retained in the",
        "developed models' (Results). Table 1 AAG on the millidose day",
        "median 0.73 g/L (IQR 0.65-1.00)."
      ),
      source_name = "AAG"
    )
  )

  compartmentData <- list(
    central = list(analyte = "docetaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "docetaxel", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    age_median = "72 years",
    age_iqr = "69-77 years",
    weight_mean = "85.5 kg",
    weight_95ci = "69.4-100.8 kg",
    bsa_mean = "2.02 m^2 (95% CI 1.79-2.23)",
    sex_female_pct = 0,
    disease_state = "Prostate cancer eligible for docetaxel 75 mg/m^2 as routine care.",
    dose_range = paste(
      "Single intravenous docetaxel millidose of 1000 ug (bolus), given 1-7",
      "days before the first therapeutic docetaxel infusion."
    ),
    regions = "Netherlands (Meander Medical Center, Amersfoort).",
    notes = paste(
      "Microdoce proof-of-concept study (EudraCT 2016-003785-77), 2017-2020;",
      "30 patients enrolled, one excluded after an allergic reaction to the",
      "millidose (Table 1 reports n = 29). Samples 0.25, 1, 2, 4, 6 and 8 h",
      "after dosing; LC-MS/MS, LLOQ 0.01 ug/L. Creatinine mean 90.9 umol/L,",
      "bilirubin mean 8.8 umol/L (Table 1). NONMEM 7.4, FOCE-I. The",
      "therapeutic-dose data of the same patients are a separate model,",
      "HeermavanVoss_2024_docetaxel_therapeutic."
    )
  )

  ini({
    # Typical values for a 70 kg patient (Table 2, 'Estimate for the millidose').
    lcl <- log(49.6); label("Clearance at 70 kg, CL (L/h)") # Table 2: Clearance 49.6 L/H (RSE 27%)
    lvc <- log(6.7); label("Central volume of distribution at 70 kg, Vc (L)") # Table 2: Central volume 6,7 L (RSE 62%); decimal comma read as 6.7
    lq <- log(39.2); label("Inter-compartmental clearance at 70 kg, Q (L/h)") # Table 2: Inter-compartmental clearance 39.2 L/H (RSE 74%)
    lvp <- log(109); label("Peripheral volume of distribution at 70 kg, Vp (L)") # Table 2: Peripheral volume 109 L (RSE 50%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL and Q (unitless)") # Methods: allometric coefficient 0.75 for flow parameters
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on Vc and Vp (unitless)") # Methods: allometric coefficient 1 for volume parameters

    # Log-normal IIV; Table 2 percentages read as CV, omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.757^2) # Table 2: IIV in clearance 75.7% (RSE 15%)
    etalvp ~ log(1 + 0.765^2) # Table 2: IIV in peripheral volume 76.5% (RSE 27%)

    propSd <- 0.370; label("Proportional residual error (fraction)") # Table 2: Residual error 37.0% (RSE 34%)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L as reported.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
