vanHoogdalem_2020_busulfan <- function() {
  description <- "Two-compartment IV PK model for busulfan in paediatric Fanconi anaemia patients undergoing haematopoietic cell transplantation, with linear fat-free-mass (FFM) scaling of all clearances and volumes (van Hoogdalem 2020)."
  reference <- "van Hoogdalem MW, Emoto C, Fukuda T, Mizuno T, Mehta PA, Vinks AA. Population pharmacokinetic modelling of busulfan and the influence of body composition in paediatric Fanconi anaemia patients. Br J Clin Pharmacol. 2020;86(5):933-943. doi:10.1111/bcp.14202"
  vignette <- "vanHoogdalem_2020_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass (Al-Sallami 2015 paediatric extension of the Janmahasatian 2005 semi-mechanistic equation, derived from total body mass, height, age and sex)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at enrolment. van Hoogdalem 2020 Equations 3 (male) and 4 (female): FFM = [a + (1 - a) / (1 + (AGE / b)^-c)] * WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT), HT in metres, with a = 0.88, b = 13.4, c = 12.7, WHSmax = 42.92, WHS50 = 30.93 for males and a = 1.11, b = 7.1, c = 1.1, WHSmax = 37.99, WHS50 = 35.98 for females. Reference FFMstd = 56.1 kg (standard adult male, 176 cm and 70 kg; Results section 3.2). The equation was developed in subjects aged 3-29 years (Discussion).",
      source_name = "FFM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body mass (TBM)",
      units = "kg",
      type = "continuous",
      notes = "Allometric TBM scaling defined the base model (exponents fixed to 0.75 / 1); FFM replaced TBM in the final model (dOFV = -28.3). TBM is an input to the FFM equation but is not referenced by the final model itself."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a covariate on the PK parameters and not retained (Results 3.2, Discussion). Used only to compute FFM."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on clearance and not retained (Discussion). Used only to select the sex-specific FFM equation."
    )
  )

  population <- list(
    n_subjects = 29,
    n_studies = 1,
    n_observations = 200,
    age_range = "4.3-22.4 years",
    age_median = "8.0 years",
    weight_range = "8.4-88.7 kg",
    weight_median = "19.9 kg",
    height_range = "86.7-178.8 cm",
    height_median = "120.6 cm",
    bmi_median = "16.0 kg/m^2",
    sex_female_pct = 55.2,
    disease_state = "Fanconi anaemia patients undergoing haematopoietic cell transplantation (busulfan-containing conditioning, phase 2 study NCT01082133); 1 overweight and 4 obese patients by the Cole age- and sex-specific BMI cut-offs",
    dose_range = "First dose of 0.6 or 0.8 mg/kg total body mass IV over 2 h (12-hour dosing interval); median initial dose 13.6 mg",
    regions = "USA (Cincinnati Children's Hospital Medical Center)",
    notes = "van Hoogdalem 2020 Table 1 and Results 3.1. Samples at 0, 0.25, 0.5, 1.5, 2, 3 and 4 h after the end of the first infusion; all but one patient younger than 18 years. Assay LLOQ 125 ng/mL."
  )

  ini({
    # Typical values standardised to FFMstd = 56.1 kg (standard adult male,
    # 70 kg and 176 cm). Table 2 labels the unit as 'per 70 kg', which is the
    # total body mass that corresponds to the FFMstd reference.
    lcl <- log(12.6);  label("Clearance at FFM = 56.1 kg (L/h)")                       # Table 2, CL = 12.6 L/h
    lvc <- log(36.4);  label("Central volume of distribution at FFM = 56.1 kg (L)")    # Table 2, V1 = 36.4 L
    lq  <- log(11.8);  label("Intercompartmental clearance at FFM = 56.1 kg (L/h)")    # Table 2, Q = 11.8 L/h
    lvp <- log(7.73);  label("Peripheral volume of distribution at FFM = 56.1 kg (L)") # Table 2, V2 = 7.73 L

    # IIV on CL and V1 (exponential model, Equation 1). CV% converted as
    # omega^2 = log(CV^2 + 1). IIV on Q and V2 was fixed to 0 in the source
    # and is therefore omitted.
    etalcl ~ 0.0259055  # Table 2, IIV for CL 16.2 CV%
    etalvc ~ 0.0138280  # Table 2, IIV for V1 11.8 CV%

    propSd <- 0.0409;  label("Proportional residual error (fraction)")                 # Table 2, proportional error 4.09 CV%
  })

  model({
    # Linear FFM scaling (Equation 2 with the exponent fixed to 1 for both
    # clearances and volumes; Results 3.2), reference FFMstd = 56.1 kg.
    size <- FFM / 56.1

    cl <- exp(lcl + etalcl) * size
    vc <- exp(lvc + etalvc) * size
    q  <- exp(lq) * size
    vp <- exp(lvp) * size

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
