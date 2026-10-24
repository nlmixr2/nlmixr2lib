Fukaya_2023_glucarpidase_healthy <- function() {
  description <- "One-compartment population PK model for glucarpidase (carboxypeptidase G2, CPG2) given as a 5-min IV infusion to healthy Japanese male volunteers (Fukaya 2023, phase 1 study), with power functions of body surface area on clearance and volume and an additive residual error."
  reference <- "Fukaya Y, Kimura T, Hamada Y, Yoshimura K, Hiraga H, Yuza Y, Ogawa A, Hara J, Koh K, Kikuta A, Koga Y, Kawamoto H. Development of a population pharmacokinetics and pharmacodynamics model of glucarpidase rescue treatment after high-dose methotrexate therapy. Front Oncol. 2023;13:1003633. doi:10.3389/fonc.2023.1003633"
  vignette <- "Fukaya_2023_glucarpidase"
  units <- list(
    time = "h",
    dosing = "mg (glucarpidase protein mass; the paper doses in U/kg and does not state the U-to-mg conversion)",
    concentration = "mg/L"
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as an un-normalised power function, tvX * BSA^theta (Fukaya 2023 Table 4). The paper does not name the BSA formula; the DuBois formula at the phase 1 mean height and weight (169.6 cm, 60.9 kg) reproduces the reported post-hoc mean CL and V.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "glucarpidase", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "20-40 years",
    age_median = "24.9 years (mean)",
    weight_range = "54.3-66.9 kg",
    weight_median = "60.9 kg (mean)",
    height_range = "162.6-184.8 cm (mean 169.6 cm)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers",
    dose_range = "Glucarpidase 20 or 50 U/kg (n = 8 each) as a 5-min IV infusion, repeated 48 h after the first dose",
    regions = "Japan (Hamamatsu University School of Medicine)",
    notes = "Phase 1 open-label randomized two-dose study (JMA-IIA00078), November 2011 to January 2012; 192 plasma CPG2 concentrations by ELISA. Demographics from Fukaya 2023 Section 3.1."
  )

  ini({
    lcl <- log(0.0590); label("Glucarpidase clearance at BSA = 1 m^2 (L/h)") # Table 4, phase 1, tvCLCPG2 = 0.0590 L/h
    lvc <- log(0.957); label("Glucarpidase volume of distribution at BSA = 1 m^2 (L)") # Table 4, phase 1, tvVCPG2 = 0.957 L
    e_bsa_cl <- 3.227; label("Power exponent of BSA on clearance (unitless)") # Table 4, phase 1, theta1 = 3.227
    e_bsa_vc <- 2.251; label("Power exponent of BSA on volume (unitless)") # Table 4, phase 1, theta2 = 2.251

    # IIV reported as CV%; omega^2 = log(CV^2 + 1)
    etalcl ~ 0.0010235 # Table 4, phase 1, CL CV 3.2%
    etalvc ~ 0.0040876 # Table 4, phase 1, V CV 6.4%

    addSd <- 0.170; label("Additive residual error (mg/L)") # Table 4, phase 1, residual variability 0.170 mg/L
  })

  model({
    cl <- exp(lcl + etalcl) * BSA^e_bsa_cl
    vc <- exp(lvc + etalvc) * BSA^e_bsa_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
