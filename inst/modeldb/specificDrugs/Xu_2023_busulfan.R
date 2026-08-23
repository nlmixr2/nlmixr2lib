Xu_2023_busulfan <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for busulfan in Chinese adult and paediatric",
    "patients undergoing allogeneic hematopoietic stem cell transplantation (Xu 2023), fitted jointly",
    "to paired plasma and saliva concentrations collected after the first 2-hour infusion. Saliva is",
    "not a separate kinetic compartment: the salivary busulfan concentration is the central-compartment",
    "(plasma) concentration multiplied by an estimated saliva:plasma scale factor of 0.88, a structure",
    "the authors selected over a kinetically distinct saliva compartment (dOFV = -82.52). Clearance and",
    "the central volume of distribution both increase with body surface area (power exponents 0.99 and",
    "1.03 referenced to the population median BSA of 1.69 m^2); the volume additionally decreases with",
    "alkaline phosphatase (power exponent -0.20 referenced to the median ALP of 74 U/L). Separate",
    "proportional residual errors apply to plasma (12.92%) and saliva (22.50%).",
    sep = " "
  )
  reference <- paste(
    "Xu B, Yang T, Zhou J, Zheng Y, Wang J, Liu Q, Li D, Zhang Y, Liu M, Wu X.",
    "Saliva as a noninvasive sampling matrix for therapeutic drug monitoring of intravenous busulfan",
    "in Chinese patients undergoing hematopoietic stem cell transplantation: A prospective population",
    "pharmacokinetic and simulation study. CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1238-1249.",
    "doi:10.1002/psp4.13004",
    sep = " "
  )
  vignette <- "Xu_2023_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    # The only ODE state. Saliva is NOT a compartment in this model: Csaliva is
    # an algebraic rescaling of the plasma concentration by fsaliva, the
    # structure the authors selected over a kinetically distinct saliva
    # compartment (Xu 2023 Results / Figure 1).
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled on both CL and Vd, referenced to the population median 1.69 m^2",
        "(Xu 2023 Table 1; range 0.92-2.14 m^2). The BSA computation formula (DuBois / Mosteller /",
        "Haycock) is not stated in the paper - record as unspecified. BSA was the single best",
        "body-size descriptor among WT, HT, BMI, IBW and AIBW in the univariate body-size screen",
        "(Xu 2023 Results, 'PopPK model')."
      ),
      source_name        = "BSA"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase activity at baseline (liver-function marker).",
      units              = "U/L (the paper prints IU/L; the two are used interchangeably)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled on Vd only, referenced to the population median 74 IU/L",
        "(Xu 2023 Table 1; range 36-253 IU/L). The estimated exponent is negative (-0.20), i.e.",
        "lower ALP gives a larger apparent volume; the authors attribute this to reduced albumin",
        "synthesis raising the free busulfan fraction (Xu 2023 Discussion)."
      ),
      source_name        = "ALP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 66,
    n_studies      = 1,
    age_range      = "6-63 years (54 adults, 12 children aged 6-17 years)",
    age_median     = "42 years",
    weight_range   = "25-90 kg",
    weight_median  = "60.70 kg",
    height_range   = "123-190 cm",
    height_median  = "160 cm",
    bsa_range      = "0.92-2.14 m^2",
    bsa_median     = "1.69 m^2",
    alp_range      = "36-253 IU/L",
    alp_median     = "74 IU/L",
    sex_female_pct = 30.30,
    race_ethnicity = c(Asian = 100),
    disease_state  = "hematologic malignancies receiving a busulfan-based myeloablative conditioning regimen prior to allogeneic hematopoietic stem cell transplantation",
    dose_range     = "0.8 mg/kg (adjusted ideal body weight) intravenously every 6 h as a 2-h infusion for 4 consecutive days",
    regions        = "China (Fujian Medical University Union Hospital, Nov 2020 - Nov 2021)",
    sampling       = "paired plasma and saliva at 0, 0.5, 2, 2.5, 3, 4 and 6 h after the start of the first dose; 406 plasma and 387 evaluable saliva samples",
    notes          = "Xu 2023 Table 1 (demographics) and Methods 'Subjects' / 'Study design' / 'Sample collection and bio-analytical assay'."
  )

  ini({
    # Structural parameters - typical values at BSA = 1.69 m^2 and ALP = 74 IU/L
    lcl <- log(8.24); label("Clearance (L/hr)")                            # Table 3 final model: theta_CL = 8.24 L/h (RSE 3%; bootstrap 7.75-8.77)
    lvc <- log(29.70); label("Central volume of distribution (L)")         # Table 3 final model: theta_Vd = 29.70 L (RSE 2%; bootstrap 28.32-31.14)

    # Saliva:plasma equilibrium concentration scale factor. The authors' final
    # structure assigns saliva to the (scaled) central compartment rather than a
    # kinetically distinct saliva compartment; the alternative structure was
    # rejected on RSE and OFV grounds (Xu 2023 Discussion).
    lfsaliva <- log(0.88); label("Saliva:plasma concentration scale factor (unitless)")  # Table 3 final model: theta_Scale-factor = 0.88 (RSE 2%; bootstrap 0.84-0.91)

    # Covariate effects - power exponents, Table 3 footnote equations
    e_bsa_cl <- 0.99; label("Power exponent for BSA on CL (unitless)")     # Table 3 final model: theta_BSA-CL = 0.99 (RSE 12%; bootstrap 0.76-1.27)
    e_bsa_vc <- 1.03; label("Power exponent for BSA on Vd (unitless)")     # Table 3 final model: theta_BSA-Vd = 1.03 (RSE 10%; bootstrap 0.77-1.22)
    e_alp_vc <- -0.20; label("Power exponent for ALP on Vd (unitless)")    # Table 3 final model: theta_ALP-Vd = -0.20 (RSE 26%; bootstrap -0.31 to -0.09)

    # IIV - exponential on CL and Vd (Xu 2023 Methods: P_i = theta_P * exp(eta_i)).
    # Table 3 reports IIV as a percentage; converted to the log-scale variance via
    # omega^2 = log(CV^2 + 1). No IIV is carried on the saliva scale factor: the
    # authors report 100% shrinkage on it and dropped it (Xu 2023 Discussion).
    etalcl ~ 0.0472764  # log(0.2200^2 + 1) = 0.0472764; Table 3 final model IIV_CL = 22.00% (RSE 8%, shrinkage 4%)
    etalvc ~ 0.0183272  # log(0.1360^2 + 1) = 0.0183272; Table 3 final model IIV_Vd = 13.60% (RSE 18%, shrinkage 14%)

    # Residual error - separate proportional models for the two matrices
    propSd <- 0.1292; label("Proportional residual SD for plasma Cc (fraction)")             # Table 3 final model: epsilon_plasma = 12.92% (RSE 17%, shrinkage 10%)
    propSd_Csaliva <- 0.2250; label("Proportional residual SD for saliva Csaliva (fraction)") # Table 3 final model: epsilon_saliva = 22.50% (RSE 22%, shrinkage 3%)
  })

  model({
    # 1. Covariate reference values (Xu 2023 Table 1 population medians, used as
    #    the normalisation constants printed in the Table 3 footnote equations)
    bsa_ref <- 1.69  # m^2
    alp_ref <- 74    # IU/L

    # 2. Individual parameters
    #    CL = theta_CL * (BSA/1.69)^theta_BSA_CL * exp(eta)
    #    Vd = theta_Vd * (BSA/1.69)^theta_BSA_Vd * (ALP/74)^theta_ALP_Vd * exp(eta)
    cl <- exp(lcl + etalcl) * (BSA / bsa_ref)^e_bsa_cl
    vc <- exp(lvc + etalvc) * (BSA / bsa_ref)^e_bsa_vc * (ALP / alp_ref)^e_alp_vc

    # Saliva:plasma scale factor (no IIV; see ini() comment)
    fsaliva <- exp(lfsaliva)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system - one compartment, IV infusion (NONMEM ADVAN1 TRANS2)
    d/dt(central) <- -kel * central

    # 5. Observations and error
    #    Plasma is the central-compartment concentration; saliva is that same
    #    concentration scaled by fsaliva (Xu 2023 Figure 1 conceptual model).
    Cc <- central / vc
    Csaliva <- fsaliva * Cc

    Cc ~ prop(propSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}
