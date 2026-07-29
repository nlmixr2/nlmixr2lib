Ji_2017_vancomycin <- function() {
  description <- "One-compartment IV (intermittent-infusion) population PK model for vancomycin in Chinese adult patients (Ji 2017). Clearance is scaled by raw Cockcroft-Gault creatinine clearance (centered linear term, reference 80 mL/min) and by age (power of (75/age), reference 75 years); the volume of distribution is a single typical value. Developed from steady-state trough therapeutic-drug-monitoring data."
  reference <- "Ji XW, Ji SM, He XR, Zhu X, Chen R, Lu W. Influences of renal function descriptors on population pharmacokinetic modeling of vancomycin in Chinese adult patients. Acta Pharmacol Sin. 2018;39(2):286-293. doi:10.1038/aps.2017.57 (published online 24 Aug 2017)"
  vignette <- "Ji_2017_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLcr. Computed by the Cockcroft-Gault equation (Ji 2017 Eq 6-7) in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Enters CL as a centered linear term referenced to 80 mL/min: TVCL is multiplied by (1 + 0.00842 * (CRCL - 80)) (Ji 2017 Eq 16). Model-building median 58.02 mL/min (range 5.45-224.0; Table 1).",
      source_name        = "CLcr"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ji 2017 Eq 16: CL is multiplied by (75/AGE)^0.08143, a power term referenced to 75 years, so typical clearance decreases as age rises above 75 (consistent with the paper's statement that vancomycin excretion declines as renal function diminishes with age). Model-building median age 78 years (range 42-95; Table 1).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 160L,
    n_studies        = 1L,
    age_range        = "42-95 years",
    age_median       = "78 years",
    weight_range     = "38-90 kg",
    weight_median    = "65 kg",
    sex_female_pct   = 33.75,
    race_ethnicity   = "Chinese adults",
    disease_state    = "Hospitalized Chinese adults receiving vancomycin for gram-positive infection; not on renal replacement therapy",
    dose_range       = "1000 mg every 12 h intravenous infusion",
    regions          = "China (Beijing Hospital, Beijing)",
    renal_function   = "Cockcroft-Gault creatinine clearance median 58.02 mL/min (range 5.45-224.0); serum creatinine median 75 umol/L (range 24-893)",
    n_concentrations = 251L,
    notes            = "Model-building cohort demographics from Ji 2017 Table 1 (n=160; a further 58 patients were held out for external validation, total n=218). Trough (Cmin) samples were collected before the next dose; each patient contributed at least one sample (median 2, range 1-17). Vancomycin concentrations measured by TDx-FLx fluorescence polarization immunoassay (LOQ 2.0 mg/L). Estimation in NONMEM 7 (FOCEI)."
  )

  ini({
    # Structural parameters (Ji 2017 Table 3 final-model column and Eq 16-17).
    # The CL reference subject is a 75-year-old with CLcr = 80 mL/min.
    lcl <- log(2.829); label("Clearance at age 75 y and CRCL 80 mL/min (L/h)") # Ji 2017 Eq 16 / Table 3: CL = 2.829 L/h
    lvc <- log(52.14); label("Volume of distribution (L)")                     # Ji 2017 Eq 17 / Table 3: V = 52.14 L

    # Covariate effects on CL (Ji 2017 Eq 16):
    #   CL = 2.829 * (1 + 0.00842 * (CLcr - 80)) * (75/Age)^0.08143 * exp(eta1)
    e_crcl_cl <- 0.00842; label("Linear CRCL slope per (CRCL - 80) mL/min on CL") # Ji 2017 Eq 16 / Table 3 theta_CLcr_CL = 0.00842
    e_age_cl  <- 0.08143; label("Power exponent on (75/AGE) for CL")             # Ji 2017 Eq 16 + Abstract + Results text = 0.08143 (Table 3 prints 0.8143 / bootstrap 0.8373 -- a dropped-leading-zero typo; see vignette Assumptions and deviations)

    # Inter-individual variability (exponential IIV, Eq 1). Ji 2017 Results text
    # reports the eta variances directly: omega^2 = 0.1051 (CL) and 0.083 (V).
    etalcl ~ 0.1051 # Ji 2017 Results: variance of eta1 on CL = 0.1051 (sqrt = 32.42%, matching Table 3 IIV_CL)
    etalvc ~ 0.083  # Ji 2017 Results: variance of eta2 on V  = 0.083  (sqrt = 28.8%, matching Table 3 IIV_V 28.87%)

    # Combined proportional + additive residual error (Eq 2; Table 3 final model).
    propSd <- 0.2679; label("Proportional residual error (fraction)") # Ji 2017 Table 3: residual CV = 26.79%
    addSd  <- 2.647;  label("Additive residual error (mg/L)")         # Ji 2017 Table 3: additive residual SD = 2.647 (table labels the column ng/mL, but the study concentration unit is mg/L throughout -- see vignette Assumptions and deviations)
  })
  model({
    # Individual clearance: typical value scaled by the CRCL centered-linear term
    # and the age power term (both on CL per Ji 2017 Eq 16), wrapped in
    # exp(etalcl) for log-normal IIV. Reference CRCL = 80 mL/min, age = 75 years.
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 80)) * (75 / AGE)^e_age_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment IV disposition. Vancomycin was given as an intermittent IV
    # infusion directly into the central compartment, so there is no absorption
    # compartment or absorption rate constant (none is estimated in Table 3).
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
