HeermavanVoss_2024_docetaxel_therapeutic <- function() {
  description <- paste(
    "Two-compartment population PK model for therapeutic intravenous",
    "docetaxel (75 mg/m^2 as a 60-minute infusion) in 29 men with prostate",
    "cancer, fitted separately from the millidose data of the same patients.",
    "Clearance and inter-compartmental clearance scale with (WT/70)^0.75 and",
    "both volumes with (WT/70)^1 (fixed exponents). Inter-individual",
    "variability on CL and the central volume; proportional residual error",
    "with a separate, larger proportional error for samples drawn during the",
    "infusion."
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
    ),
    TINF = list(
      description = "Duration of the current (most recent) intravenous infusion, in hours.",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not a covariate on a PK parameter: it selects the residual error.",
        "Some therapeutic-dose samples were drawn during the infusion by",
        "mistake, and the source estimated a separate proportional error",
        "for them (Results; Table 2 'Residual error for samples drawn during",
        "infusion' 98.8%). rxode2 does not expose the infusion duration to",
        "model(), so a sample is flagged as during the infusion when",
        "tad(central) < TINF. The study infusion lasted 60 min (TINF = 1)."
      ),
      source_name = "(not named; samples drawn during infusion)"
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
        "developed models' (Results). Table 1 AAG on the therapeutic-dose",
        "day median 0.72 g/L (IQR 0.63-0.88)."
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
    disease_state = "Prostate cancer treated with docetaxel 75 mg/m^2 as routine care.",
    dose_range = paste(
      "Docetaxel 75 mg/m^2 as a 60-minute intravenous infusion; Table 1",
      "mean dose 151.2 mg (95% CI 133.6-168.0)."
    ),
    regions = "Netherlands (Meander Medical Center, Amersfoort).",
    notes = paste(
      "Microdoce proof-of-concept study (EudraCT 2016-003785-77), 2017-2020;",
      "30 patients enrolled, one excluded after an allergic reaction to the",
      "preceding millidose (Table 1 reports n = 29). Samples 0.25, 1, 2, 4, 6",
      "and 8 h after the end of the infusion, some drawn during the infusion",
      "by mistake; LC-MS/MS, LLOQ 0.01 ug/L. NONMEM 7.4, FOCE-I. The",
      "millidose data of the same patients are a separate model,",
      "HeermavanVoss_2024_docetaxel_millidose."
    )
  )

  ini({
    # Typical values for a 70 kg patient (Table 2, 'Estimate for the therapeutic dose').
    lcl <- log(34.6); label("Clearance at 70 kg, CL (L/h)") # Table 2: Clearance 34.6 L/H (RSE 34%)
    lvc <- log(8.00); label("Central volume of distribution at 70 kg, Vc (L)") # Table 2: Central volume 8.00 L (RSE 39%)
    lq <- log(4.63); label("Inter-compartmental clearance at 70 kg, Q (L/h)") # Table 2: Inter-compartmental clearance 4.63 L/H (RSE 20%)
    lvp <- log(22.2); label("Peripheral volume of distribution at 70 kg, Vp (L)") # Table 2: Peripheral volume 22.2 L (RSE 32%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL and Q (unitless)") # Methods: allometric coefficient 0.75 for flow parameters
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on Vc and Vp (unitless)") # Methods: allometric coefficient 1 for volume parameters

    # Log-normal IIV; Table 2 percentages read as CV, omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.617^2) # Table 2: IIV in clearance 61.7% (RSE 61%)
    etalvc ~ log(1 + 0.692^2) # Table 2: IIV in central volume 69.2% (RSE 87%)

    propSd <- 0.471; label("Proportional residual error, samples after the end of infusion (fraction)") # Table 2: Residual error 47.1% (RSE 21%)
    propSd_infusion <- 0.988; label("Proportional residual error, samples drawn during the infusion (fraction)") # Table 2: Residual error for samples drawn during infusion 98.8% (RSE 57%)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L as reported.
    Cc <- central / vc * 1000

    # tad(central) runs from the start of the infusion, so a sample is during
    # the infusion while it is shorter than the infusion duration.
    infusing <- 0
    if (tad(central) < TINF) infusing <- 1
    propSdEff <- propSd * (1 - infusing) + propSd_infusion * infusing
    Cc ~ prop(propSdEff)
  })
}
