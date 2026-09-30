Mileva_2021_doxycycline_rabbit <- function() {
  description <- "Preclinical (rabbit). One-compartment population PK model with first-order absorption and elimination for doxycycline hyclate given as a single 5 mg/kg oral dose in hard gelatin capsules to mature (5-month-old) and immature (70-day-old) New Zealand x Californian broiler rabbits, parameterised per kg body weight with no covariates (Mileva 2021)"
  reference <- "Mileva R, Rusenov A, Milanova A. Population Pharmacokinetic Modelling of Orally Administered Doxycycline to Rabbits at Different Ages. Antibiotics (Basel). 2021;10(3):310. doi:10.3390/antibiotics10030310"
  vignette <- "Mileva_2021_doxycycline_rabbit"
  # Every structural parameter is normalised to body weight exactly as
  # published: the dose is given in mg/kg, compartment amounts are carried in
  # mg/kg, V/F in L/kg and CL/F in L/kg/h, so central/vc is mg/L = ug/mL and no
  # unit conversion is applied.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "doxycycline", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a covariate but not retained in the final model (Mileva 2021 Results section 2 and Methods section 4.5: 'did not improve the tested model'). Mean +/- SD 3.55 +/- 0.30 kg (mature) and 2.03 +/- 0.21 kg (immature), Methods section 4.2. Because V/F and CL/F are expressed per kg, whole-animal V/F and CL/F are proportional to body weight.",
      source_name = "body weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a covariate (age group: 70-day-old immature vs 5-month-old mature rabbits, i.e. about 0.19 vs 0.42 years) but not retained in the final model (Mileva 2021 Results section 2).",
      source_name = "age"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose biochemistry (Mileva 2021 Supplementary Table S2: 65.02 +/- 5.7 g/L immature, 69.98 +/- 4.81 g/L mature). Tested but not retained in the final model.",
      source_name = "Total protein"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose biochemistry (Mileva 2021 Supplementary Table S2: 28.42 +/- 1.11 g/L immature, 31.4 +/- 2.65 g/L mature). Tested but not retained in the final model.",
      source_name = "Albumin"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose biochemistry in IU (Mileva 2021 Supplementary Table S2: 73.25 +/- 11.17 immature, 34.17 +/- 11.29 mature). Tested but not retained in the final model.",
      source_name = "ALT (UI)"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose biochemistry in IU (Mileva 2021 Supplementary Table S2: 44.75 +/- 7.69 immature, 34.67 +/- 14.54 mature). Tested but not retained in the final model.",
      source_name = "AST (UI)"
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose biochemistry in IU (Mileva 2021 Supplementary Table S2: 294.17 +/- 70.99 immature, 319.67 +/- 42.87 mature). Tested but not retained in the final model.",
      source_name = "LDH (UI)"
    )
  )

  population <- list(
    species = "rabbit (crossbred New Zealand x Californian broiler)",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "70 days (immature, n = 10) and 5 months (mature, n = 6)",
    weight_range = "mean +/- SD 2.03 +/- 0.21 kg (immature) and 3.55 +/- 0.30 kg (mature)",
    sex_female_pct = NA_real_,
    disease_state = "clinically healthy",
    dose_range = "5 mg/kg doxycycline hyclate, single oral dose in a hard gelatin capsule",
    regions = "Bulgaria (Trakia University, Stara Zagora)",
    sampling = "Mature: 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 24 h. Immature: two alternating sub-groups sampled at 0.5, 2, 4, 8, 12 h or at 1, 3, 6, 10, 24 h. Free (TFA-precipitated) plasma doxycycline by HPLC-PDA, LOQ 0.15 ug/mL; censored samples were excluded from the fit.",
    notes = "18 rabbits enrolled (6 mature, 12 immature; Methods section 4.2); two immature rabbits that expelled the capsule were excluded, leaving 16 in the analysis (Figure 1 caption, Supplementary Table S1). Sex was not reported. Fit in Phoenix NLME 8.3 with FOCE-ELS."
  )

  ini({
    lka <- log(0.257); label("Absorption rate constant ka (1/h)") # Table 1, tvka = 0.257 1/h
    lcl <- log(1.473); label("Apparent clearance CL/F per kg body weight (L/kg/h)") # Table 1, tVCl = 1.473 L/kg/h
    lvc <- log(4.429); label("Apparent volume of distribution V/F per kg body weight (L/kg)") # Table 1, tvV = 4.429 L/kg

    etalka ~ 0.103 # Table 2, eta ka variance 0.103 (BSV 32.96% CV)
    etalcl ~ 0.033 # Table 2, eta Cl variance 0.033 (BSV 18.18% CV)
    etalvc ~ 0.392 # Table 2, eta V variance 0.392 (BSV 69.30% CV)

    propSd <- 0.368; label("Proportional residual error (fraction)") # Table 1, stdev0 = 0.368; multiplicative error per Equation 8
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
