Steffens_2025_amikacin <- function() {
  description <- "One-compartment IV population PK model for amikacin in Brazilian hospitalized adult and pediatric patients undergoing therapeutic drug monitoring, with an exponential creatinine-clearance effect on CL (Steffens 2025)"
  reference <- "Steffens NA, Zimmermann ES, Azeredo FJ, Linden R, Finatto LJ, Hahn RZ, Schwarzbold AV, Pacheco LS, Brucker N. Therapeutic Drug Monitoring-Based Population Pharmacokinetics of Amikacin in Patients at a Teaching Hospital. Antibiotics (Basel). 2025;14(6):531. doi:10.3390/antibiotics14060531"
  vignette <- "Steffens_2025_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Estimated by the Cockcroft-Gault equation in raw mL/min (Steffens 2025 Section 4.3), NOT BSA-normalized to mL/min/1.73 m^2. Stored under the canonical CRCL column per inst/references/covariate-columns.md, which admits raw mL/min when the source paper applies no BSA normalization; the same convention is used by the sibling amikacin model Delattre_2010_amikacin.R. Entered UNCENTERED on the log scale: log(CL) = log(1.49) + 0.004 * CRCL, so exp(lcl) is CL extrapolated to CRCL = 0 rather than CL at a typical patient. See the ini() comment on e_crcl_cl for the reconciliation that establishes the uncentered form. Population median 79.01 mL/min, range 12.97-517.97 (Steffens 2025 Table 1).",
      source_name        = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a continuous covariate by forward-inclusion / backward-elimination (Steffens 2025 Section 4.4.1) but NOT retained in the final model; only CRCL improved the fit (Steffens 2025 Section 2.2). Median 69.40 kg, range 15.60-143.80."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained (Steffens 2025 Section 2.2). Median 51 years, range 4-75. The authors note likely collinearity between age and creatinine clearance (Steffens 2025 Discussion)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Steffens 2025 Section 2.2). Median 23.20 kg/m^2, range 10.60-50.70. The authors note likely collinearity between weight and BMI (Steffens 2025 Discussion)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as a categorical covariate but not retained (Steffens 2025 Sections 2.2 and 4.4.1). 6 of 39 patients (15.4%) were female."
    ),
    DIALYSIS = list(
      description = "Dialysis status indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as a categorical covariate but not retained (Steffens 2025 Section 4.4.1). No per-patient counts are reported."
    )
  )

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "4-75 years",
    age_median     = "51 years",
    weight_range   = "15.60-143.80 kg",
    weight_median  = "69.40 kg",
    sex_female_pct = 15.4,
    race_ethnicity = "Not reported (Brazilian teaching-hospital population)",
    disease_state  = "Critically and non-critically ill hospitalized patients treated with amikacin for at least three days for Gram-negative and multidrug-resistant infections; most commonly Klebsiella pneumoniae (n = 20) and Pseudomonas aeruginosa (n = 7)",
    dose_range     = "225-1500 mg amikacin IV, administered every 12, 24, 48 or 72 h; most common regimen 1000 mg q24h (n = 15, 35.7%)",
    regions        = "Brazil (Hospital Universitario de Santa Maria, Rio Grande do Sul)",
    renal_function = "Cockcroft-Gault creatinine clearance median 79.01 mL/min (range 12.97-517.97), raw mL/min and not BSA-normalized; 7 patients below 30 mL/min and 10 above 120 mL/min",
    notes          = "Baseline demographics per Steffens 2025 Table 1. Single-center prospective observational study, May 2018 - February 2020; 43 patients enrolled and 4 excluded for inappropriate collection time or missing data. Cohort is predominantly adult (28 adults, 9 elderly) but includes 2 pediatric patients, so the age range extends to 4 years. Burns, pregnancy and refusal of consent were exclusion criteria. A total of 113 amikacin concentrations (53 peak, 60 trough; 2-6 samples per subject) were collected by non-routine therapeutic drug monitoring at steady state after at least three days of therapy: trough 30 min before a dose and peak 30 min after the end of infusion. Mean +/- SD concentrations were 41.96 +/- 20.20 ug/mL (peak) and 8.75 +/- 15.38 ug/mL (trough). Fit in Monolix 2024R1 by SAEM; final estimates confirmed by a 500-sample nonparametric bootstrap."
  )

  ini({
    # Structural parameters, Steffens 2025 Table 2 "Final Estimate" column.
    #
    # NOTE ON THE CL INTERCEPT: exp(lcl) is CL extrapolated to CRCL = 0, NOT CL
    # at a typical patient. The paper prints the covariate coefficient
    # ("beta CrCl on Cl" = 0.004) but never writes the covariate equation, so
    # the centering had to be established from the paper's own numbers. Monolix
    # enters a continuous covariate UNCENTERED by default and names the
    # coefficient exactly as printed here (beta_Cl_CrCl); under that reading the
    # final model reproduces the paper's own covariate-free base model at the
    # population median CRCL:
    #   1.49 * exp(0.004 * 79.01) = 2.044 L/h  vs  base-model CL = 2.04 L/h
    #     (Steffens 2025 Section 2.2, median CRCL from Table 1)
    # A median-centered reading would instead put CL = 1.49 L/h at the median,
    # a 27% unexplained drop from the base model. See the vignette Errata.
    lcl <- log(1.49);  label("Clearance extrapolated to CRCL = 0 (L/h)")     # Steffens 2025 Table 2: Cl = 1.49 L/h (RSE 12.96%); intercept of log(CL) = log(Cl) + 0.004 * CRCL
    lvc <- log(23.18); label("Central volume of distribution (L)")           # Steffens 2025 Table 2: Vd = 23.18 L (RSE 23.73%)

    # Covariate effect: exponential (log-linear), uncentered, on CL.
    e_crcl_cl <- 0.004; label("Exponential CRCL effect on CL (per mL/min)")  # Steffens 2025 Table 2: beta CrCl on Cl = 0.004 (RSE 57.15%)

    # Inter-individual variability. Table 2 prints omega as the log-scale SD
    # with the lognormal %CV in parentheses; ini() takes the VARIANCE, so each
    # entry below is omega^2. Confirmed against the printed %CV:
    #   sqrt(exp(0.67^2) - 1) = 75.3%  vs printed 74.8%
    #   sqrt(exp(0.47^2) - 1) = 49.7%  vs printed 49.89%
    # (the small residual is the 2-decimal rounding of omega itself; the
    # variance reading would give 97.7% and 77.5%, which the paper does not print).
    etalcl ~ 0.4489 # 0.67^2; Steffens 2025 Table 2: omega Cl = 0.67 (74.8% CV), RSE 17.01%
    etalvc ~ 0.2209 # 0.47^2; Steffens 2025 Table 2: omega Vd = 0.47 (49.89% CV), RSE 37.69%

    # Residual error. Section 2.2 records that a combined additive-plus-
    # proportional error model was TESTED, but the final model reported in
    # Table 2 and in the Abstract carries a proportional term only.
    propSd <- 0.38; label("Proportional residual error (fraction)") # Steffens 2025 Table 2: proportional error model = 0.38 (RSE 14.86%)
  })
  model({
    # Individual PK parameters. CRCL enters CL exponentially and uncentered,
    # so CL rises from 1.49 L/h at anuria to 2.04 L/h at the population median
    # CRCL of 79.01 mL/min.
    cl <- exp(lcl + e_crcl_cl * CRCL + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc is mg/L, which is numerically
    # identical to the ug/mL used throughout Steffens 2025.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
