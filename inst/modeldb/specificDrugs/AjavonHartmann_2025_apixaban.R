AjavonHartmann_2025_apixaban <- function() {
  description <- "Population PK model for oral apixaban in adults and pediatric patients, including children with congenital or acquired heart disease (SAXOPHONE) (Ajavon-Hartmann 2025): two-compartment with first-order absorption and elimination, dose-dependent relative bioavailability, body-weight power effects on CL/F and Vc/F, a fixed literature CYP maturation function on CL/F, age on Q/F, a 9-month-to-18-year age-group effect and an ALL/LL patient-type effect on ka, and a congenital-or-acquired-heart-disease effect on CL/F, with a zero-intercept linear mixed-effects PK/PD relationship to anti-factor Xa activity."
  reference <- "Ajavon-Hartmann A, Jarugula P, Back H, Obianom O, Ludwig E, Crevar C, Marchisin D, Wang Z, Chen W, He B, Payne RM, Burns KM, Glatz AC, Male C, Monagle P, Perera V, Murthy B, Merali S. Characterization of Apixaban Pharmacokinetics/Pharmacodynamics and Dose Assessment in Pediatric Patients with Congenital or Acquired Heart Disease. Clin Pharmacol Ther. 2025;118:365-377. doi:10.1002/cpt.3689"
  vignette <- "AjavonHartmann_2025_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (exponent 0.68) and Vc/F (exponent 0.882), normalised to the 67.4-kg reference patient of the apixaban pediatric analyses (Table 1 equations). Time-fixed at baseline (source column BBWT).",
      source_name = "BBWT"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Three uses: (1) power effect (AGE / 25)^0.751 on Q/F, 25 years being the reference age; (2) the fixed literature CYP maturation function FBirth + (1 - FBirth) * AGE^0.83 / (AGE^0.83 + 0.244^0.83) on CL/F; (3) derivation of the age-group indicator PED9M (1 for 9 months <= AGE < 18 years, else 0) that scales ka. Postnatal age in years (Table 1 footnote: 'AGEi is the postnatal age in years').",
      source_name = "AGE"
    ),
    DIS_ALL = list(
      description = "Pediatric patient with acute lymphoblastic leukemia or lymphoblastic lymphoma treated with asparaginase (PREVAPIX-ALL patient type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy adults, pediatric patients at risk of venous or arterial thromboembolism from phase I studies CV185079 and CV185118, and pediatric patients with congenital or acquired heart disease from SAXOPHONE)",
      notes = "Exponential effect on ka (exp(-0.988)); unlike the Jarugula 2025 fit, backward elimination removed the patient-type effect on Vc/F (Figure S1). Labelled PTYPE_ALL; Table 1 footnote: 'Pediatric patient with ALL at risk of VTE denotes pediatric patients with ALL or lymphoblastic lymphoma at risk of VTE treated with asparaginase', so the indicator pools leukemia and lymphoblastic lymphoma and is 1 for every PREVAPIX-ALL (CV185155) patient.",
      source_name = "PTYPE"
    ),
    DIS_HEART_CONGENITAL_ACQUIRED = list(
      description = "Pediatric patient with congenital or acquired heart disease requiring anticoagulation for thromboembolism prevention (SAXOPHONE patient type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy adults, pediatric patients at risk of venous or arterial thromboembolism, and pediatric patients with ALL or lymphoblastic lymphoma)",
      notes = "Exponential effect on CL/F (exp(-0.185), i.e. 16.9% lower CL/F at the same weight; Table 1 theta17). Labelled PTYPE_CAHD; 1 for every SAXOPHONE (CV185362) patient. The authors attribute the lower clearance to most patients having had a Fontan procedure (Discussion).",
      source_name = "PTYPE"
    ),
    DOSE_APIXABAN_MG = list(
      description = "Administered apixaban dose amount",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent relative bioavailability: F1 = 1 for doses <= 2.5 mg and 1 - expit(-0.457) * ((DOSE - 2.5) / 47.5)^0.89 above it (Table 1 equations; DOSE is 'the dose regimen assigned by weight tier'). Unlike the ALL/LL regimen, the CAHD regimen reaches 3, 4 and 5 mg in the heavier weight tiers, so the term is active there (F1 about 0.99, 0.98 and 0.97). Must agree with the event-table amt on the dose records.",
      source_name = "DOSE"
    ),
    FORM_CAPSULE = list(
      description = "Sprinkle-capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet or oral solution, relative bioavailability 1)",
      notes = "1 = the apixaban sprinkle capsule. Multiplies relative bioavailability by the fixed value 1.1 (Table 1 theta14 '1.1 FIXED'; FFORM equals theta14 for the sprinkle capsule and 1 for others).",
      source_name = "FFORM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 602,
    n_studies = 12,
    age_range = "Pediatric 0.044-18 years (Table S1); SAXOPHONE 0.378-18 years; adult ages not tabulated",
    age_median = "Pediatric 6.73 years (Table S1, overall pediatric); SAXOPHONE 7.5 years",
    weight_range = "Pediatric 3-133 kg (Table S1); SAXOPHONE 6.1-133 kg; adult weights not tabulated",
    weight_median = "Pediatric 22.6 kg (Table S1, overall pediatric); SAXOPHONE 23.6 kg",
    sex_female_pct = 47.3,
    race_ethnicity = c(
      White = 78.5,
      Black = 6.44,
      Asian = 7.43,
      `American Indian or Alaska Native` = 0.248,
      `Native Hawaiian or Pacific Islander` = 0.248,
      Other = 7.18
    ),
    disease_state = "Pooled healthy adults (198, eight adult studies), pediatric patients at risk of venous or arterial thromboembolism (56, phase I studies CV185079 and CV185118), pediatric patients with acute lymphoblastic leukemia or lymphoblastic lymphoma receiving asparaginase (224, PREVAPIX-ALL CV185155), and pediatric patients aged 28 days to <18 years with congenital or acquired heart disease requiring chronic anticoagulation for thromboembolism prevention (124, phase II SAXOPHONE CV185362, NCT02981472).",
    dose_range = "SAXOPHONE fixed-dose weight-tiered oral twice-daily regimen: 1, 1.5, 2, 3, 4 and 5 mg for 6 to <9, 9 to <12, 12 to <18, 18 to <25, 25 to <35 and >=35 kg (implied by Table 2 dose / CL/F; proposed regimen in the Results); adult studies up to 50 mg.",
    regions = "Multinational",
    n_observations = 6098,
    notes = "Sex and race percentages are for the 404 pediatric patients of Table S1; adult demographics are not reported. Of the 6098 concentrations, 447 were from SAXOPHONE. The PK/PD (anti-factor Xa) layer was fitted to 390 AXA values from 123 PREVAPIX-ALL and SAXOPHONE patients (Results; the Methods quote 541 concentrations from 253 patients for the PK/PD dataset before exclusions) measured with a drug-calibrated assay in ng/mL."
  )

  ini({
    # Structural parameters: Table 1 (final PPK model), 'Estimate' column.
    lka <- log(0.525)
    label("Absorption rate constant for adults (ka, 1/h)") # Table 1 theta1 = 0.525
    lcl <- log(4.7)
    label("Apparent clearance for a 67.4-kg fully matured reference patient (CL/F, L/h)") # Table 1 theta2 = 4.7
    lvc <- log(32.5)
    label("Apparent central volume for a 67.4-kg reference patient (Vc/F, L)") # Table 1 theta3 = 32.5
    lq <- log(1.43)
    label("Apparent intercompartmental clearance at the 25-year reference age (Q/F, L/h)") # Table 1 theta4 = 1.43
    lvp <- log(20.8)
    label("Apparent peripheral volume of distribution (Vp/F, L)") # Table 1 theta5 = 20.8

    # Dose-dependent relative bioavailability (Table 1 equations).
    lgamma <- log(0.89)
    label("Shape factor for the reduction in F1 at doses above 2.5 mg (theta6, unitless)") # Table 1 theta6 = 0.89
    logitifdepot50 <- -0.457
    label("Logit of the maximum reduction in F1 (reached at a 50-mg dose) (theta7, unitless)") # Table 1 theta7 = -0.457
    e_form_capsule_fdepot <- fixed(1.1)
    label("Relative bioavailability of the sprinkle-capsule formulation (theta14, ratio)") # Table 1 theta14 = 1.1 FIXED

    # Covariate effects: Table 1.
    e_wt_vc <- 0.882
    label("Power exponent on (WT / 67.4) for Vc/F (theta8, unitless)") # Table 1 theta8 = 0.882
    e_wt_cl <- 0.68
    label("Power exponent on (WT / 67.4) for CL/F (theta9, unitless)") # Table 1 theta9 = 0.68
    e_age_q <- 0.751
    label("Power exponent on (AGE / 25) for Q/F (theta10, unitless)") # Table 1 theta10 = 0.751
    e_ped9m_ka <- 1.23
    label("Fractional change in ka for patients aged 9 months to <18 years (theta11, unitless)") # Table 1 theta11 = 1.23
    e_dis_all_ka <- -0.988
    label("Exponential effect of the ALL/LL patient type on ka (theta16, unitless)") # Table 1 theta16 = -0.988
    e_dis_heart_congenital_acquired_cl <- -0.185
    label("Exponential effect of the congenital-or-acquired-heart-disease patient type on CL/F (theta17, unitless)") # Table 1 theta17 = -0.185

    # Non-renal CL/F maturation function taken from the literature (Table 1
    # thetas 12, 13 and 15, all fixed).
    lhill <- fixed(log(0.83))
    label("Shape (Hill) parameter of the non-renal CL/F ontogeny function (theta12, unitless)") # Table 1 theta12 = 0.83 FIXED
    lpna50 <- fixed(log(0.244))
    label("Age at half-maximal non-renal CL/F maturation, AGE50 (theta13, years)") # Table 1 theta13 = 0.244 FIXED
    fbirth_cl <- fixed(0.05)
    label("Fraction of adult non-renal CL/F activity at birth, FBirth (theta15, fraction)") # Table 1 theta15 = 0.05 FIXED

    # Linear mixed-effects PK/PD model, Table 3.
    slope_antiFXa <- 0.916
    label("Slope of anti-factor Xa activity vs apixaban plasma concentration (unitless, ng/mL per ng/mL)") # Table 3 Slope = 0.916

    # IIV: Table 1, variances (SD in parentheses in the table). Footnote
    # correlations 0.809 (ZCL:ZVC) and 0.912 (ZQ:ZVP).
    etalka ~ 0.399 # Table 1 'ZKA' omega 1,1 = 0.399
    etalcl + etalvc ~ c(0.149, 0.0838, 0.0719) # Table 1 'ZCL' 0.149, 'ZVC:ZCL' 0.0838, 'ZVC' 0.0719
    etalvp + etalq ~ c(1.39, 1.07, 0.989) # Table 1 'ZVP' 1.39, 'ZQ:ZVP' 1.07, 'ZQ' 0.989

    # PK/PD random slope: Table 3 'Random effects, Standard deviation' = 0.0528,
    # additive on the slope (linear mixed-effects regression); 0.0528^2.
    etaslope_antiFXa ~ 0.00278784

    # Residual error: Table 1, additive on log concentration, variance 0.108
    # (SD 0.328).
    expSd <- 0.328
    label("Additive residual SD on log-transformed apixaban concentration") # Table 1 sigma 1,1 = 0.108 (SD 0.328)

    addSd_antiFXa <- 12.8
    label("Additive residual SD on anti-factor Xa activity (ng/mL)") # Table 3 'Residual' = 12.8
  })
  model({
    # Age-group indicator PED9M: 1 for 9 months to <18 years (Table 1 footnote).
    ped9m <- (AGE >= 0.75) * (AGE < 18)

    # Absorption rate constant, Table 1 equations.
    ka <- exp(lka + etalka) * (1 + e_ped9m_ka * ped9m) * exp(e_dis_all_ka * DIS_ALL)

    # Non-renal maturation of CL/F (maturation function MF, Table 1 footnote).
    hill <- exp(lhill)
    pna50 <- exp(lpna50)
    fmat <- fbirth_cl + (1 - fbirth_cl) * AGE^hill / (AGE^hill + pna50^hill)

    cl <- exp(lcl + etalcl) * (WT / 67.4)^e_wt_cl * fmat *
      exp(e_dis_heart_congenital_acquired_cl * DIS_HEART_CONGENITAL_ACQUIRED)
    vc <- exp(lvc + etalvc) * (WT / 67.4)^e_wt_vc
    q <- exp(lq + etalq) * (AGE / 25)^e_age_q
    vp <- exp(lvp + etalvp)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # Dose-dependent relative bioavailability anchored at 1 for doses <= 2.5 mg;
    # the excess is clamped at zero so lower doses cannot raise a negative
    # number to a fractional power.
    gamma <- exp(lgamma)
    ifdepot50 <- expit(logitifdepot50)
    doseExcess <- (DOSE_APIXABAN_MG > 2.5) * (DOSE_APIXABAN_MG - 2.5)
    frel <- 1 - ifdepot50 * (doseExcess / 47.5)^gamma
    fform <- 1 + (e_form_capsule_fdepot - 1) * FORM_CAPSULE
    f(depot) <- frel * fform

    # Dose in mg and volume in L give mg/L; x1000 reports ng/mL.
    Cc <- 1000 * central / vc

    # Zero-intercept linear mixed-effects PK/PD (Discussion: 'intercept fixed
    # to zero and IIV in the slope parameter'); drug-calibrated activity in ng/mL.
    slope_i <- slope_antiFXa + etaslope_antiFXa
    antiFXa <- slope_i * Cc

    Cc ~ lnorm(expSd)
    antiFXa ~ add(addSd_antiFXa)
  })
}
