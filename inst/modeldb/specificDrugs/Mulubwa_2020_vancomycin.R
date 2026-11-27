Mulubwa_2020_vancomycin <- function() {
  description <- "One-compartment population PK model for intravenous vancomycin in preterm neonates (Mulubwa 2020). Clearance and volume scale allometrically with current body weight (fixed exponents 0.75 and 1, reference 1.48 kg); log-normal between-subject variability on CL and V and exponential residual error. Parameters were estimated by maximum a posteriori Bayesian estimation in Monolix from trough concentrations only."
  reference <- "Mulubwa M, Griesel HA, Mugabo P, Dippenaar R, van Wyk L. Assessment of Vancomycin Pharmacokinetics and Dose Regimen Optimisation in Preterm Neonates. Drugs R D. 2020;20(2):105-113. doi:10.1007/s40268-020-00302-7"
  vignette <- "Mulubwa_2020_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin was infused intravenously over 60 min (Methods, 'Vancomycin
  # Dosing, Blood Sampling and Plasma Assay'), so the dose enters `central`.
  # Samples were collected in SST tubes, allowed to clot and centrifuged
  # (serum), although the paper calls the result a 'plasma concentration'.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight (at the time of the study)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 1.48 kg, the cohort median current weight (Table 1; Eqs. 2-3; Table 2 footnote a). Enters CL with fixed exponent 0.75 and V with fixed exponent 1. Current weight was preferred over BSA as the body-size descriptor (Results, 'Pharmacokinetic Model').",
      source_name = "Weight"
    )
  )

  # Screened one at a time in the covariate analysis but without significant
  # effect on CL or V (Results, 'Pharmacokinetic Model'); documentation only.
  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and V; not significant (Results). Cohort median 33 (30-34.7) weeks (Table 1).",
      source_name = "PMA"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not significant (Results). Cohort median 31 (23-34) weeks (Table 1).",
      source_name = "gestational age"
    ),
    PNA = list(
      description = "Postnatal (current) age",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not significant (Results). Cohort median 14 (3-58) days (Table 1).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Tested; not significant (Results). 11 of 19 female (Table 1).",
      source_name = "sex"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the allometric size descriptor in place of current weight; current weight described body size better (Results). Cohort median 0.131 (0.091-0.18) m^2 (Table 1).",
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the 17 neonates with values; not significant (Results). Median 27 (19-37) g/L (Table 1).",
      source_name = "albumin"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the 16 neonates with values; not significant (Results). Median 51 (26-74) umol/L (Table 1).",
      source_name = "serum creatinine"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, eGFR = 0.33 x height / serum creatinine (Eq. 1)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the 16 neonates with values; not significant (Results). Table 1 prints a median of 0.28 (0.12-0.51) with units mL/min/1.73 m^2; Eq. 1 with height in cm and creatinine in umol/L gives values in those magnitudes, i.e. the printed value is not on the conventional mL/min/1.73 m^2 scale.",
      source_name = "eGFR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    n_observations = 45L,
    age_range = "3-58 days postnatal",
    age_median = "14 days postnatal",
    weight_range = "0.925-2.62 kg (current weight)",
    weight_median = "1.48 kg (current weight)",
    sex_female_pct = 57.9,
    ga_range = "23-34 weeks (median 31)",
    pma_range = "30-34.7 weeks (median 33)",
    birth_weight = "median 1.27 kg (0.63-2.69)",
    disease_state = "Preterm neonates in neonatal intensive care with suspected nosocomial (MRSA) infection; neonates with renal dysfunction or chromosomal abnormalities excluded. Seven had respiratory distress syndrome, two jaundice and one necrotising enterocolitis.",
    dose_range = "15 mg/kg every 12 h as a 60-min intravenous infusion (Aspen-Vancomycin or Sandoz-Vancocin CP, randomised)",
    regions = "South Africa (Netcare Blaauwberg and N1 City hospitals, Cape Town)",
    notes = "Randomised double-blind brand-comparison study (Methods; Table 1). Trough samples only, drawn before the 3rd, 6th and 9th doses (2-3 per neonate). Estimation by MAP Bayesian estimation (SAEM) in Monolix 2018R2 with priors from a published neonatal vancomycin model fitted to peak and trough data from 54 neonates (reference 18 of the paper)."
  )

  ini({
    # Structural parameters for a typical neonate of 1.48 kg current weight
    # (Table 2 footnote a; Eqs. 2-3).
    lcl <- log(0.102); label("Clearance at WT = 1.48 kg (L/h)") # Table 2 CL 0.102 L/h (RSE 8.19%); Eq. 2
    lvc <- log(0.884); label("Volume of distribution at WT = 1.48 kg (L)") # Table 2 Vd 0.884 L (RSE 24%); Eq. 3

    # Allometric exponents held constant (Methods 'Pharmacokinetic Analysis';
    # Results: estimating them did not improve precision).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/1.48) on CL (unitless)") # Methods; Eq. 2 exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/1.48) on V (unitless)") # Methods; Eq. 3 exponent 1

    # Between-subject variability (Table 2, as %CV; footnote b defines
    # CV = sqrt(exp(SD^2) - 1), so omega^2 = log(CV^2 + 1)).
    etalcl ~ 0.0015987 # Table 2 BSV CL 4 %CV (RSE 63.4%): log(1 + 0.04^2)
    etalvc ~ 0.25588 # Table 2 BSV Vd 54 %CV (RSE 21.4%): log(1 + 0.54^2)

    # Monolix exponential error model: log(y) = log(f) + a * e.
    expSd <- 0.284; label("Exponential (log-scale) residual error SD (unitless)") # Table 2 residual error 'Exponential' 0.284 (RSE 13.2%)
  })
  model({
    # Individual parameters (Eqs. 2-3)
    cl <- exp(lcl + etalcl) * (WT / 1.48)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 1.48)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L, so central / vc is mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
