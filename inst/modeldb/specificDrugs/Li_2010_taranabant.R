Li_2010_taranabant <- function() {
  description <- "Three-compartment population PK model for oral taranabant in healthy and obese adults (Li 2010)"
  reference <- "Li X, Nielsen J, Cirincione B, Li H, Addy C, Wagner J, Hartford A, Erondu N, Gantz I, Morgan J, Stone J. Development of a Population Pharmacokinetic Model for Taranabant, a Cannabinoid-1 Receptor Inverse Agonist. AAPS J. 2010;12(4):537-547. doi:10.1208/s12248-010-9212-2"
  vignette <- "Li_2010_taranabant"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in Phase 2 obese cohort (Li 2010 Methods). Power-form effect on CL/F and on V4/F with reference 31.5 kg/m^2.",
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on Q4/F and on V4/F (additive age contribution) with reference age 39 years.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Additive 643 L shift on V4/F for female subjects. Source column SEXF i in Li 2010 Equation 3 reads 1 for females.",
      source_name        = "SEXF"
    ),
    CRCL = list(
      description        = "Creatinine clearance, raw (not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in Phase 2 obese cohort (Li 2010 Methods). Linear-deviation effect on CL/F (slope 0.0668 L/h per mL/min) and on V4/F (slope 12.5 L per mL/min) centered at 80.6 mL/min. Cockcroft-Gault is presumed but not stated in the paper.",
      source_name        = "CrCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 572,
    n_studies      = 13,
    age_range      = "18-69 years",
    age_median     = "39 years",
    weight_range   = "50-151 kg",
    weight_median  = "88.9 kg",
    sex_female_pct = 64.3,
    race_ethnicity = c(White = 61.9, Black = 17.0, Hispanic = 14.3, Asian = 6.1, Other = 0.7),
    disease_state  = "Pooled Phase 1 healthy volunteers (n = 187) and Phase 2 obese subjects with BMI 30-43 kg/m^2 (n = 385); 12 Phase 1 studies and 1 Phase 2 study",
    dose_range     = "0.5-8 mg oral (single and multiple daily doses; multiple doses >= 10 mg excluded from analysis)",
    regions        = "Multinational including a Japanese Phase 1 sub-study",
    renal_function = "Includes mild and moderate renal insufficiency from Phase 1 study 029 (CRCL range 5.1-171 mL/min). Subjects with moderate hepatic insufficiency were excluded.",
    notes          = "Baseline demographics from Li 2010 Table II (pooled column). 6,834 plasma concentrations contributed to the fit."
  )

  ini({
    # Structural PK parameters (Li 2010 Table III, final estimates)
    # Reference covariate values embedded in the typical-value equations:
    # BMI_ref = 31.5 kg/m^2, AGE_ref = 39 years, CRCL_ref = 80.6 mL/min.

    lka  <- log(2.73);  label("Absorption rate constant (1/h)")                          # Table III; ka = 2.73 1/h
    lcl  <- log(25.4);  label("Apparent clearance base coefficient for BMI power term (L/h)")  # Table III Eq. 1; CL/F = 25.4 L/h at reference covariates
    lvc  <- log(175);   label("Apparent central volume of distribution (L)")             # Table III; V2/F = 175 L
    lvp  <- log(273);   label("Apparent first peripheral volume of distribution (L)")    # Table III; V3/F = 273 L (no IIV estimated)
    lvp2 <- log(2130);  label("Apparent second peripheral volume base coefficient for BMI power term (L)")  # Table III Eq. 3; V4/F base = 2130 L
    lq   <- log(24.6);  label("Apparent intercompartmental clearance to peripheral1 (L/h)")  # Table III; Q3/F = 24.6 L/h
    lq2  <- log(30.4);  label("Apparent intercompartmental clearance base coefficient for Q4 age power term (L/h)")  # Table III Eq. 2; Q4/F base = 30.4 L/h at AGE_ref

    # Paper-specific additive age contribution to V4/F (intercept 752 L scaled by (AGE/39)^e_age_vp2)
    lvp2_age <- log(752); label("Additive age-intercept contribution to V4/F (L)")      # Table III Eq. 3; intercept for age on V4/F = 752 L

    # Covariate effects (Li 2010 Table III, final estimates)
    e_bmi_cl  <- -1.11;   label("BMI power exponent on CL/F (unitless)")                # Table III; reference BMI 31.5 kg/m^2
    e_crcl_cl <-  0.0668; label("Linear CrCL slope on CL/F (L/h per mL/min)")           # Table III; centered at 80.6 mL/min
    e_bmi_vp2 <-  1.38;   label("BMI power exponent on V4/F base term (unitless)")      # Table III; reference BMI 31.5 kg/m^2
    e_age_vp2 <-  2.10;   label("AGE power exponent on V4/F age term (unitless)")       # Table III; reference age 39 years
    e_sexf_vp2 <- 643;    label("Additive shift on V4/F for female subjects (L)")       # Table III; SEXF = 1 for females
    e_crcl_vp2 <- 12.5;   label("Linear CrCL slope on V4/F (L per mL/min)")             # Table III; centered at 80.6 mL/min
    e_age_q2  <-  0.357;  label("AGE power exponent on Q4/F (unitless)")                # Table III; reference age 39 years

    # IIV: omega^2 = log(CV^2 + 1) where CV is reported in Table III as %CV.
    # ka (73.0%), Q3 (37.2%), Q4 (35.9%) share a 3 x 3 correlated block;
    # CL (44.1%), Vc (29.7%), Vp2 (33.3%) are diagonal. V3 IIV not estimated (NE).
    etalka + etalq + etalq2 ~ c(0.4272,
                                0.0901, 0.1297,
                                0.0676, 0.0764, 0.1212)  # Table III diagonals + covariances
    etalcl  ~ 0.1779  # 44.1% CV
    etalvc  ~ 0.0846  # 29.7% CV
    etalvp2 ~ 0.1051  # 33.3% CV

    # Residual error (Li 2010 Table III, Phase 2 prop+add model).
    # Phase 1 used a TSLD-stratified proportional-only model:
    #   prop variance 0.407 for TSLD < 2 h, 0.0526 for TSLD >= 2 h (no additive).
    # Phase 2 (clinical cohort) used a single additive + proportional model
    # without time stratification because few samples were collected before 2 h.
    # The model file encodes the Phase 2 structure as a representative
    # cross-population error. propSd is unit-free; additive variance
    # 0.00385 was reported in nM (per LLOQ of 0.1 nM giving 71.5% CV).
    # The additive component is dropped because the model concentration
    # unit is ug/mL (= mg/L) rather than nM; see vignette Errata.
    propSd <- 0.355; label("Proportional residual error (Phase 2 fraction; sqrt(0.126))")  # Table III; Phase 2 prop variance 0.126
  })
  model({
    # Typical clearance with BMI power effect + linear CrCL deviation (Li 2010 Eq. 1)
    tvcl <- exp(lcl) * (BMI / 31.5)^e_bmi_cl + e_crcl_cl * (CRCL - 80.6)
    cl   <- tvcl * exp(etalcl)

    # Typical V4/F (second peripheral) is an ADDITIVE sum of four covariate
    # contributions, not a multiplicative scaling (Li 2010 Eq. 3):
    #   V4 = 2130*(BMI/31.5)^1.38 + 752*(AGE/39)^2.10 + 643*SEXF + 12.5*(CRCL-80.6)
    tvvp2 <- exp(lvp2)     * (BMI / 31.5)^e_bmi_vp2 +
             exp(lvp2_age) * (AGE / 39  )^e_age_vp2 +
             e_sexf_vp2 * SEXF +
             e_crcl_vp2 * (CRCL - 80.6)
    vp2  <- tvvp2 * exp(etalvp2)

    # Typical Q4/F with AGE power effect (Li 2010 Eq. 2)
    tvq2 <- exp(lq2) * (AGE / 39)^e_age_q2
    q2   <- tvq2 * exp(etalq2)

    # Remaining structural parameters (no covariates; ka/Q3 share IIV block with Q4)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)            # V3/F: no IIV estimated (Table III)
    q  <- exp(lq  + etalq)

    # Micro-constants for explicit ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Observation: dose in mg, vc in L gives mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
