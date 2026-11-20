Bae_2026_busulfan <- function() {
  description <- "One-compartment population PK model for intravenous busulfan in pediatric and young adult allogeneic hematopoietic stem cell transplant recipients (Bae 2026), with elimination scaled by an enzyme-turnover auto-inhibition compartment: busulfan concentration inhibits the zero-order synthesis of a normalised enzyme pool (IC50 3250 ug/L, turnover rate 0.0404 1/h) so that clearance falls over a 4-day once-daily course. Volume scales with BSA as a power function (exponent 1.21) and clearance linearly with BSA, both normalised to 1.73 m^2, with an exponential age-maturation of clearance reaching 50% at 0.374 years."
  reference <- paste(
    "Bae S, Choi JY, Hong KT, Kim BK, Park HJ, Yu KS, Rhee SJ, Kang HJ.",
    "Clinical Impact of Maximum Plasma Concentration of Busulfan in Pediatric",
    "Patients Undergoing Allogenic Hematopoietic Stem Cell Transplantation.",
    "Clin Pharmacol Ther. 2026;120(3):702-712. doi:10.1002/cpt.70337.",
    sep = " "
  )
  vignette <- "Bae_2026_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # The enzyme state is dimensionless: it is normalised to 1 at baseline
  # (Methods: 'The initial amount of enzyme is set to 1').
  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE),
    enzyme = list(
      analyte = "busulfan-metabolising enzyme pool (interpreted by the authors as glutathione S-transferase)",
      units = "(fraction of baseline)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalised to 1.73 m^2 on both V (power, exponent 1.21) and CL (linear). The paper does not state the BSA formula used.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the exponential maturation of clearance, CL * (1 - exp(-0.693 / thalf_age_cl * AGE)); the maturation factor is 0.5 at AGE = 0.374 years and approaches 1 in older children and adults.",
      source_name = "Age"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 334,
    n_studies = 1,
    n_observations = 5357,
    age_range = "0.6-30.3 years",
    age_mean = "8.8 years (3-h infusion group) and 10.0 years (6-h infusion group)",
    weight_range = "7.4-91.3 kg",
    weight_mean = "32.3 kg (3-h group) and 33.5 kg (6-h group)",
    bsa_range = "0.26-2.13 m^2",
    sex_female_pct = 40.7,
    race_ethnicity = c(Asian_Korean = 100),
    disease_state = "Pediatric and young adult allogeneic hematopoietic stem cell transplant recipients (high-risk malignancy 75.7%, immunodeficiency 10.5%, other rare genetic disease 6.9%, congenital hematologic disease 4.2%, histiocytosis 2.7%) receiving busulfan-based myeloablative conditioning.",
    dose_range = "Once-daily intravenous busulfan for 4 days, infused over 3 h (2009-2013, n = 122) or 6 h (2014-2020, n = 212), with therapeutic drug monitoring targeting a cumulative 4-day AUC of 75,000 ug*h/L; mean total dose 416-427 mg/m^2.",
    regions = "Korea (single centre: Seoul National University Children's Hospital)",
    notes = "Bae 2026 Table 1 and Methods. Retrospective TDM data: samples at the end of infusion and 1, 2 and 4 h after, sometimes also before the next day's dose; 5,357 concentrations by LC-MS/MS. NONMEM 7.5.1. Only 30 patients (9.0%) were aged 18 years or older. The site is inferred from the institution; race is not tabulated."
  )

  ini({
    # Volume: V = theta1 * (BSA / 1.73)^theta2 (Table 3 footnote a).
    lvc <- log(42)
    label("Volume of distribution at BSA 1.73 m^2 (L)") # Table 3 'theta1; typical V value with BSA 1.73 m2' = 42 (RSE 2%)
    e_bsa_vc <- 1.21
    label("Power exponent of BSA/1.73 on V (unitless)") # Table 3 'theta2; BSA exponent for V' = 1.21 (RSE 3%)

    # Clearance: CL = theta3 * (BSA / 1.73) * (1 - exp(-0.693 / theta4 * Age)) (Table 3 footnote b).
    lcl <- log(11)
    label("Clearance of a fully mature subject at BSA 1.73 m^2 and baseline enzyme (L/h)") # Table 3 'theta3; typical CL value with BSA 1.73 m2' = 11 (RSE 2%)
    thalf_age_cl <- 0.374
    label("Age at 50% of mature clearance in an exponential maturation (years)") # Table 3 'theta4; Age reaching 50% of adult maturation' = 0.374 (RSE 13%)

    # Enzyme turnover auto-inhibition (Eqs 1-2 and Figure S2).
    lic50 <- log(3250)
    label("Busulfan concentration giving 50% inhibition of enzyme synthesis, IC50 (ug/L)") # Table 3 'IC50; Half maximal inhibitory concentration' = 3,250 ug/L (RSE 8%)
    lkenz <- log(0.0404)
    label("Enzyme turnover rate constant, K_ENZ (1/h)") # Table 3 'K ENZ; Enzyme turnover rate' = 0.0404 1/h (RSE 17%)

    # IIV reported as CV%: omega^2 = log(1 + CV^2). The off-diagonal is read
    # as a covariance on the omega scale (Table 3 labels it a correlation);
    # see the vignette Assumptions and deviations section.
    etalcl + etalvc ~ c(0.04769, 0.0226, 0.02915) # Table 3 'IIV for CL (%)' = 22.1 -> log(1 + 0.221^2) = 0.04769; 'Correlation between IIV on CL and V' = 0.0226; 'IIV for V (%)' = 17.2 -> log(1 + 0.172^2) = 0.02915

    propSd <- 0.103
    label("Proportional residual error (fraction)") # Table 3 'Proportional residual error' = 0.103 (RSE 3%)
  })

  model({
    # Age maturation of clearance (Table 3 footnote b, constant 0.693 as printed).
    mat_cl <- 1 - exp(-0.693 / thalf_age_cl * AGE)

    cl <- exp(lcl + etalcl) * (BSA / 1.73) * mat_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.73)^e_bsa_vc
    ic50 <- exp(lic50)
    kenz <- exp(lkenz)

    # Dose in mg and V in L give mg/L; x 1000 gives ug/L, the unit of IC50.
    Cc <- 1000 * central / vc

    # Eq 1: elimination is first order in busulfan and proportional to the
    # relative enzyme amount.
    d/dt(central) <- -cl / vc * central * enzyme
    # Eq 2 with first-order enzyme loss kenz * enzyme (Figure S2 draws the
    # loss arrow as K_ENZ out of the enzyme compartment; the typeset Eq 2
    # omits the '* A2'). Synthesis is inhibited by busulfan concentration.
    d/dt(enzyme) <- kenz * (1 - Cc / (ic50 + Cc)) - kenz * enzyme
    enzyme(0) <- 1

    Cc ~ prop(propSd)
  })
}
