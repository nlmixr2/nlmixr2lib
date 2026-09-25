Smit_2020_vancomycin <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous vancomycin in",
    "morbidly obese adults undergoing bariatric surgery and nonobese healthy",
    "volunteers (Smit 2020). Clearance scales with total body weight by an",
    "estimated power exponent; the first peripheral volume (V2) scales",
    "linearly with total body weight; central (V1) and first peripheral",
    "volume share a linear age effect centred on 36.5 years. Between-subject",
    "variability on clearance is estimated separately for the obese and",
    "nonobese cohorts, selected by DIS_OBESE_MORBID."
  )
  reference <- paste(
    "Smit C, Wasmann RE, Goulooze SC, Wiezer MJ, van Dongen EPA, Mouton JW,",
    "Bruggemann RJM, Knibbe CAJ. Population pharmacokinetics of vancomycin",
    "in obesity: Finding the optimal dose for (morbidly) obese individuals.",
    "Br J Clin Pharmacol. 2020;86(2):303-317. doi:10.1111/bcp.14144."
  )
  vignette <- "Smit_2020_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed total body weight (TBW). Power covariate on CL",
        "(CL70kg * (TBW/70)^0.535) and linear-proportional covariate on the",
        "first peripheral volume V2 (V2_70kg;36.5yr * TBW/70), both",
        "referenced to 70 kg (Table 2)."
      ),
      source_name = "TBW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear covariate (1 + theta2 * (AGE - 36.5)) on V1 and V2 with a",
        "single shared slope theta2 = 0.0136 per year (Table 2)."
      ),
      source_name = "Age"
    ),
    DIS_OBESE_MORBID = list(
      description = "Morbidly obese cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonobese healthy volunteer)",
      notes = paste(
        "1 = morbidly obese patient undergoing bariatric surgery (BMI >= 40",
        "kg/m^2, or >= 35 kg/m^2 with comorbidities), 0 = nonobese healthy",
        "volunteer (BMI 18-25 kg/m^2). Selects the cohort-specific",
        "between-subject variability on clearance (Table 2: 'CL nonobese'",
        "5.28% FIX, 'CL obese' 24.7%). It has no effect on any typical value."
      ),
      source_name = "obese / nonobese study group (Table 2 IIV split)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 28,
    n_studies = 1,
    age_range = "20-55 years",
    age_median = "obese 38.0 years; nonobese 25.5 years",
    weight_range = "60.0-234.6 kg",
    weight_median = "obese 139.0 kg; nonobese 69.5 kg",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = paste(
      "20 morbidly obese adults (BMI 40.8-65.7 kg/m^2) undergoing laparoscopic",
      "sleeve gastrectomy or gastric bypass, and 8 nonobese healthy",
      "volunteers (BMI 20.4-25.0 kg/m^2). All had eGFR >= 60 mL/min/1.73 m^2",
      "(measured 24-h creatinine clearance median 141.4 mL/min obese, 117.9",
      "mL/min nonobese)."
    ),
    dose_range = paste(
      "Single IV infusion at 10 mg/min: obese 12.5 mg/kg TBW (maximum 2500",
      "mg), nonobese 1000 mg."
    ),
    regions = "The Netherlands (St. Antonius Hospital, Nieuwegein)",
    notes = paste(
      "Demographics from Smit 2020 Table 1 (sex distribution not reported).",
      "326 plasma samples over 48 h (11-13 per subject); 24 samples (7%)",
      "below the 1.5 mg/L limit of detection were handled with the M3",
      "method. NONMEM 7.4, FOCE-I with LAPLACIAN. Externally validated on",
      "previously published Blouin et al. data (6 obese, 4 nonobese; ref. 21)."
    )
  )

  ini({
    # Structural parameters: Smit 2020 Table 2, final model column.
    lcl <- log(5.72); label("Clearance for a 70 kg individual (L/h)") # Table 2: CL70kg = 5.72 L/h (RSE 5.0%)
    lvc <- log(16.7); label("Central volume V1 at age 36.5 years (L)") # Table 2: V1_36.5yr = 16.7 L (RSE 18%)
    lq <- log(15.8); label("Intercompartmental clearance V1-V2 (L/h)") # Table 2: Q V1-V2 = 15.8 L/h (RSE 23%)
    lvp <- log(6.98); label("First peripheral volume V2 at 70 kg and 36.5 years (L)") # Table 2: V2_70kg;36.5yr = 6.98 L (RSE 17%)
    lq2 <- log(5.21); label("Intercompartmental clearance V1-V3 (L/h)") # Table 2: Q V1-V3 = 5.21 L/h (RSE 21%)
    lvp2 <- log(19.5); label("Second peripheral volume V3 (L)") # Table 2: V3 = 19.5 L (RSE 13%)

    # Covariate effects
    e_wt_cl <- 0.535; label("Power exponent of (WT/70) on CL (unitless)") # Table 2: theta1 = 0.535 (RSE 20%)
    e_age_vc_vp <- 0.0136; label("Linear slope of (AGE - 36.5) on V1 and V2 (1/year)") # Table 2: theta2 = 0.0136 (RSE 31%), shared by V1 and V2

    # IIV (Table 2 footnote b: CV% = sqrt(exp(omega^2) - 1), so omega^2 = log(1 + CV^2)).
    # Separate CL variances per cohort; multiplexed by DIS_OBESE_MORBID in model().
    etalcl ~ fixed(0.002784) # Table 2: IIV CL nonobese = 5.28% FIX; log(1 + 0.0528^2) = 0.002784
    etalcl_obese ~ 0.05920 # Table 2: IIV CL obese = 24.7% (RSE 19%); log(1 + 0.247^2) = 0.05920
    etalvc ~ 0.18664 # Table 2: IIV V1 = 45.3% (RSE 24%); log(1 + 0.453^2) = 0.18664

    # Residual error: Table 2 footnote c, proportional error shown as sigma (SD).
    propSd <- 0.0392; label("Proportional residual error (fraction)") # Table 2: proportional error = 0.0392 (RSE 21%)
    addSd <- 1.07; label("Additive residual error (mg/L)") # Table 2: additive error = 1.07 mg/L (RSE 5.0%)
  })

  model({
    # Cohort-specific CL variability (Table 2 'CL nonobese' / 'CL obese').
    etacl_cohort <- etalcl * (1 - DIS_OBESE_MORBID) + etalcl_obese * DIS_OBESE_MORBID

    # Individual parameters (Table 2 final-model equations)
    cl <- exp(lcl + etacl_cohort) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (1 + e_age_vc_vp * (AGE - 36.5))
    q <- exp(lq)
    vp <- exp(lvp) * (WT / 70) * (1 + e_age_vc_vp * (AGE - 36.5))
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
