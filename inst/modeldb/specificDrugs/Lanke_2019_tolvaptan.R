Lanke_2019_tolvaptan <- function() {
  description <- "One-compartment population PK model for oral tolvaptan in adults with autosomal dominant polycystic kidney disease (Lanke 2019): zero-order absorption directly into the central compartment (duration D1), first-order elimination, and a dose-dependent relative bioavailability F1 = Dosemax / (Dose50 + Dose) that falls by about 37 percent from 15 mg to 120 mg per administration. eGFR (power), body weight (power) and CYP3A inhibitor coadministration (fractional change) act on CL/F; body weight (power) acts on Vd/F. Correlated IIV on CL/F, Vd/F and D1, independent IIV on Dose50 and Dosemax, and combined additive plus proportional residual error."
  reference <- paste(
    "Lanke S, Shoaf SE. (2019).",
    "Population Pharmacokinetic Analyses and Model Validation of Tolvaptan",
    "in Subjects With Autosomal Dominant Polycystic Kidney Disease.",
    "J Clin Pharmacol 59(5):763-770.",
    "doi:10.1002/jcph.1370"
  )
  vignette <- "Lanke_2019_tolvaptan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "tolvaptan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (exponent -0.53) and Vd/F (exponent 0.23), both normalised to the cohort median of 79.5 kg (Lanke 2019 Equation 1 and Table 1).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Creatinine-based CKD-EPI eGFR (Lanke 2019 Methods, Covariate Model Building). Time-varying in the source dataset; enters CL/F as a power term (exponent 0.52) normalised to the cohort median of 69.4 mL/min/1.73 m^2 (Lanke 2019 Equation 2 and Table 1).",
      source_name = "eGFR"
    ),
    CONMED_CYP3A4_INH = list(
      description = "CYP3A inhibitor coadministration indicator (1 = coadministered on this record)",
      units = "(binary)",
      type = "binary",
      reference_category = "No CYP3A inhibitor coadministered",
      notes = "Pools strong, moderate and weak CYP3A inhibitors (Lanke 2019 Methods, 131 subjects with at least one instance). Changes record by record within a subject (time-varying). Enters CL/F as theta_typ * (1 + e_cyp3a4_inh_cl * CONMED_CYP3A4_INH) with e_cyp3a4_inh_cl = -0.23 (Lanke 2019 Equation 3 and Table 3). 3.46 percent of observations carried the indicator (Table 2).",
      source_name = "Inh"
    ),
    DOSE_TOLVAPTAN_MG = list(
      description = "Tolvaptan amount of the current administration (mg)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent relative bioavailability f(central) = fdepot_dosemax / (fdepot_dose50 + DOSE_TOLVAPTAN_MG) (Lanke 2019 Equation 4). Per-administration amount, not the daily total: for a 90/30 mg split regimen the column is 90 on the morning dose record and 30 on the afternoon dose record. Place this column AFTER the event columns (id, time, evid, amt, cmt) in a data frame passed to rxSolve(), because rxode2 can read a dose-named column that precedes amt as the dose amount.",
      source_name = "Dose"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tabulated as a candidate covariate (Lanke 2019 Table 1, median 42 years, range 18-62) but not retained in the final model.",
      source_name = "Age"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tabulated as a candidate covariate (Lanke 2019 Table 1, median 25.9, range 17.31-90.8) but not retained.",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "Male",
      notes = "Tabulated as a candidate covariate (Lanke 2019 Table 2, 526 male / 565 female) but not retained.",
      source_name = "Sex"
    ),
    REGION_JAPAN = list(
      description = "Enrolled at a site in Japan",
      units = "(binary)",
      type = "binary",
      reference_category = "Non-Japan site",
      notes = "Tabulated as a candidate covariate (Lanke 2019 Table 2, 18 of 1091 subjects enrolled in Japan) but not retained.",
      source_name = "Site"
    ),
    CONMED_CYP3A4_IND = list(
      description = "CYP3A4 inducer coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "No CYP3A4 inducer coadministered",
      notes = "Excluded from the forward covariate search because fewer than 1 percent of observations (64 of 7355) followed CYP3A4 inducer coadministration (Lanke 2019 Table 2 and Discussion).",
      source_name = "Ind"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1091L,
    n_observations = 7335L,
    n_studies = 6L,
    age_range = "18-62 years",
    age_median = "42 years",
    weight_range = "42.3-164.5 kg",
    weight_median = "79.5 kg",
    sex_female_pct = 51.8,
    race_ethnicity = "Of the 1073 subjects enrolled outside Japan: 1024 White, 14 Black, 1 Hispanic, 9 Asian, 25 other; 18 subjects enrolled in Japan",
    disease_state = "Autosomal dominant polycystic kidney disease (ADPKD)",
    renal_function = "CKD-EPI eGFR median 69.4 mL/min/1.73 m^2 (range 13.62-137.37)",
    dose_range = "Oral tolvaptan tablets: single doses of 15, 30, 60 or 120 mg, 30 mg once daily, or daily split-dose regimens 15/15, 30/15, 30/30, 45/15, 60/30 and 90/30 mg (morning dose about 8 AM, second dose 8-9 h later)",
    regions = "International (1073 subjects outside Japan, 18 in Japan)",
    co_medication = "131 subjects with at least one record of strong, moderate or weak CYP3A inhibitor coadministration; 22 subjects with a CYP3A4 inducer",
    notes = "Model development data: 5 phase 2 trials and the phase 3b open-label extension NCT01214421 (Lanke 2019 Methods); dense:sparse data ratio 1:5. External validation (not used for estimation): 678 subjects, 3012 observations from the phase 3b randomized-withdrawal trial NCT02160145 at 60/30 or 90/30 mg. Food intake was not recorded in the phase 3b trials. Assay range 5.00-1000 ng/mL; concentrations below the LLOQ (1.2 percent) were excluded. Demographics from Lanke 2019 Tables 1 and 2."
  )

  ini({
    # Structural parameters -- Lanke 2019 Table 3, Final Model column.
    lcl <- log(12.6) ; label("Apparent clearance CL/F at median eGFR and WT, no CYP3A inhibitor (L/h)") # Table 3 'CL/F (L/h)' final = 12.6
    lvc <- log(110)  ; label("Apparent volume of distribution Vd/F at median WT (L)")                   # Table 3 'Vd/F (L)' final = 110
    ld1 <- log(0.58) ; label("Zero-order absorption duration D1 (h)")                                   # Table 3 'D1 (h)' final = 0.58

    # Dose-dependent relative bioavailability F1 = Dosemax / (Dose50 + Dose)
    # (Lanke 2019 Equation 4); F1 is relative to a 15 mg administration,
    # where 182 / (166 + 15) = 1.006.
    lfdepot_dose50  <- log(166) ; label("Dose-scale constant in the F1 denominator, Dose50 (mg)") # Table 3 'Dose50 (mg)' final = 166 (Results prose prints 160; see vignette)
    lfdepot_dosemax <- log(182) ; label("Dose-scale constant in the F1 numerator, Dosemax (mg)")  # Table 3 'Dosemax (mg)' final = 182

    # Covariate effects -- Lanke 2019 Table 3 final model.
    e_crcl_cl <- 0.52        ; label("Power exponent of eGFR / 69.4 on CL/F (unitless)")    # Table 3 'Power coefficient eGFR on CL/F' = 0.52
    e_wt_cl <- -0.53         ; label("Power exponent of WT / 79.5 on CL/F (unitless)")      # Table 3 'Power coefficient WT on CL/F' = -0.53
    e_cyp3a4_inh_cl <- -0.23 ; label("Fractional change in CL/F with CYP3A inhibitor (unitless)") # Table 3 'Fractional change of CL/F with Inh' = -0.23
    e_wt_vc <- 0.23          ; label("Power exponent of WT / 79.5 on Vd/F (unitless)")      # Table 3 'Power coefficient WT on Vd/F' = 0.23

    # IIV -- Table 3 'Interindividual variability' block, final model. The
    # entries are omega variances: sqrt(0.49) = 70 and sqrt(5.67) = 238
    # reproduce the Abstract's 70 and 238 percent IIV in Vd/F and D1.
    # Block order (CL, Vd, D1): var CL 0.40; cov CL-Vd -0.25; var Vd 0.49;
    # cov CL-D1 1.13; cov Vd-D1 -1.51; var D1 5.67.
    etalcl + etalvc + etald1 ~ c(
      0.40,
      -0.25, 0.49,
      1.13, -1.51, 5.67
    )
    etalfdepot_dose50 ~ 0.086  # Table 3 IIV 'Dose50' final = 0.086
    etalfdepot_dosemax ~ 0.019 # Table 3 IIV 'Dosemax' final = 0.019

    # Residual error -- Table 3 'Residual variability' final: Additive 5.3,
    # Proportional 0.18, read as NONMEM SIGMA variances (operator ruling on
    # the scale, sidecar 2026-09-29), so SD = sqrt(variance).
    addSd <- 2.3022  ; label("Additive residual error SD (ng/mL)")    # sqrt(5.3); Table 3 'Additive' final = 5.3
    propSd <- 0.4243 ; label("Proportional residual error SD (fraction)") # sqrt(0.18); Table 3 'Proportional' final = 0.18
  })

  model({
    # Covariates normalised to the cohort medians (Lanke 2019 Table 1):
    # WT 79.5 kg, eGFR 69.4 mL/min/1.73 m^2. Equations 1-3.
    cl <- exp(lcl + etalcl) *
      (CRCL / 69.4)^e_crcl_cl *
      (WT / 79.5)^e_wt_cl *
      (1 + e_cyp3a4_inh_cl * CONMED_CYP3A4_INH)
    vc <- exp(lvc + etalvc) * (WT / 79.5)^e_wt_vc
    d1 <- exp(ld1 + etald1)
    fdepot_dose50 <- exp(lfdepot_dose50 + etalfdepot_dose50)
    fdepot_dosemax <- exp(lfdepot_dosemax + etalfdepot_dosemax)

    kel <- cl / vc

    # One compartment; the dose enters central at a constant rate over d1
    # (Lanke 2019 Equation 5). Dose records need rate = -2 for dur() to apply.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # Relative bioavailability of the current administration (Equation 4).
    f(central) <- fdepot_dosemax / (fdepot_dose50 + DOSE_TOLVAPTAN_MG)

    # mg / L -> ng/mL
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
