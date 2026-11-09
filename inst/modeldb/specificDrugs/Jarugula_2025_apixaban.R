Jarugula_2025_apixaban <- function() {
  description <- "Population PK model for oral apixaban in adults and pediatric patients, including children with acute lymphoblastic leukemia or lymphoblastic lymphoma treated with asparaginase (PREVAPIX) (Jarugula 2025): two-compartment with first-order absorption and elimination, dose-dependent relative bioavailability, body-weight power effects on CL/F and Vc/F, a fixed literature CYP maturation function on CL/F, age on Q/F, a 9-month-to-18-year age-group effect on ka and patient-type (ALL/LL) effects on ka and Vc/F, with a zero-intercept linear PK/PD relationship to anti-factor Xa activity."
  reference <- "Jarugula P, Ajavon-Hartmann T, Obianom O, Ludwig E, Crevar C, Marchisin D, Wang Z, Chen W, He B, Mitchell LG, Rodriguez V, Lew G, Perera V, Murthy B, Merali S. Modeling and Simulation to Inform Apixaban Dosing in Pediatrics With Acute Lymphoblastic Leukemia or Lymphoblastic Lymphoma Treated With Asparaginase. Clin Pharmacol Ther. 2025;118(2):355-364. doi:10.1002/cpt.3616"
  vignette <- "Jarugula_2025_apixaban"
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
      notes = "Power effects on CL/F (exponent 0.669) and Vc/F (exponent 0.840), normalised to the 67.4-kg reference patient of the apixaban pediatric analyses (Table 2 equations). Time-fixed at baseline (control-stream column BBWT).",
      source_name = "BBWT"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Three uses: (1) power effect (AGE / 25)^0.719 on Q/F, 25 years being the reference age of the apixaban pediatric analyses; (2) the fixed literature CYP maturation function FBirth + (1 - FBirth) * AGE^0.83 / (AGE^0.83 + 0.244^0.83) on CL/F; (3) derivation of the age-group indicator PED9M (1 for 9 months <= AGE < 18 years, else 0) that scales ka. Postnatal age in years; the source control stream reads AGE in years throughout.",
      source_name = "AGE"
    ),
    DIS_ALL = list(
      description = "Pediatric patient with acute lymphoblastic leukemia or lymphoblastic lymphoma treated with asparaginase during induction chemotherapy (PREVAPIX patient type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adults from the eight adult studies, and pediatric patients at risk of venous or arterial thromboembolism from phase I studies CV185079 and CV185118)",
      notes = "Exponential effects on ka (exp(-1.67)) and Vc/F (exp(-0.564)). The indicator is the control stream's PTYPE_CAN ('IF(PTYPE.EQ.2) PTYPE_CAN = 1'), labelled PTYPE_ALL in Table 2; it is 1 for every PREVAPIX patient (224 of 224, Table 1) and 0 for everyone else, so it pools leukemia and lymphoblastic lymphoma and is confounded with study. The authors note the effect may partly reflect between-study and sampling-design differences rather than disease alone (Discussion).",
      source_name = "PTYPE"
    ),
    DOSE_APIXABAN_MG = list(
      description = "Administered apixaban dose amount",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent relative bioavailability: F1 = 1 for doses <= 2.5 mg and 1 - expit(-0.461) * ((DOSE - 2.5) / 47.5)^0.889 above it (Table 2 equations; control stream 'IF(DOSA.LE.2.5) THEN F1=1 ELSE ... F1=(1-IMF*((DOSA-2.5)/47.5)**GAM)'). Every pediatric weight-tier dose (0.1 to 2.5 mg) sits at or below the 2.5-mg anchor, so F1 is 1 for them. Must agree with the event-table amt on the dose records.",
      source_name = "DOSA"
    ),
    FORM_CAPSULE = list(
      description = "Sprinkle-capsule (new neonatal) formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet or oral solution, relative bioavailability 1)",
      notes = "1 = the apixaban sprinkle capsule developed for neonates and infants. Multiplies relative bioavailability by the fixed literature value 1.1 (Table 2 theta14 '1.1 FIXED'; control stream 'FFORM = 1; IF(GROUPNUM.EQ.10) FFORM=THETA(14)'). The capsule is the arm carrying the non-unit relative bioavailability; tablet and oral solution are the reference.",
      source_name = "GROUPNUM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 478,
    n_studies = 11,
    age_range = "Pediatric 0.044-18 years (Table 1); adult ages not tabulated",
    age_median = "Pediatric 6.29 years (Table 1, overall pediatric)",
    weight_range = "Pediatric 3-121 kg (Table 1); adult weights not tabulated",
    weight_median = "Pediatric 22.6 kg (Table 1, overall pediatric)",
    sex_female_pct = 45.4,
    race_ethnicity = c(
      White = 75.7,
      Black = 6.79,
      Asian = 8.57,
      `Native Hawaiian or Pacific Islander` = 0.357,
      Other = 8.57
    ),
    disease_state = "Pooled adults (198, eight adult studies), pediatric patients at risk of venous or arterial thromboembolism (56, phase I studies CV185079 and CV185118), and pediatric patients aged 1 to <18 years with newly diagnosed acute lymphoblastic leukemia or lymphoblastic lymphoma with a central venous access device receiving asparaginase during induction chemotherapy (224, phase III PREVAPIX, NCT02369653).",
    dose_range = "PREVAPIX fixed-dose weight-tiered oral regimen: 0.5, 1, 1.5, 2 and 2.5 mg twice daily for 6 to <10.5, 10.5 to <18, 18 to <25, 25 to <35 and >=35 kg; adult studies up to 50 mg.",
    regions = "Multinational",
    n_observations = 5651,
    notes = "Sex and race percentages are for the 280 pediatric patients of Table 1 only; adult demographics are not reported in this paper. Of the 5651 concentrations, 431 were from PREVAPIX. The PK/PD (anti-factor Xa) layer was fitted to 151 activity values from 130 PREVAPIX patients measured with the Diagnostica Stago anti-Xa apixaban assay (drug-calibrated, ng/mL)."
  )

  ini({
    # Structural parameters: Table 2 (final PPK model), 'Estimate' column.
    lka <- log(0.527)
    label("Absorption rate constant for adults (ka, 1/h)") # Table 2 theta1 = 0.527
    lcl <- log(4.70)
    label("Apparent clearance for a 67.4-kg fully matured reference patient (CL/F, L/h)") # Table 2 theta2 = 4.70
    lvc <- log(32.7)
    label("Apparent central volume for a 67.4-kg reference patient (Vc/F, L)") # Table 2 theta3 = 32.7
    lq <- log(1.45)
    label("Apparent intercompartmental clearance at the 25-year reference age (Q/F, L/h)") # Table 2 theta4 = 1.45
    lvp <- log(21.0)
    label("Apparent peripheral volume of distribution (Vp/F, L)") # Table 2 theta5 = 21.0

    # Dose-dependent relative bioavailability (Table 2 equations; Data S1 $PK).
    lgamma <- log(0.889)
    label("Shape factor for the reduction in F1 at doses above 2.5 mg (theta6, unitless)") # Table 2 theta6 = 0.889
    logitifdepot50 <- -0.461
    label("Logit of the maximum reduction in F1 (reached at a 50-mg dose) (theta7, unitless)") # Table 2 theta7 = -0.461
    e_form_capsule_fdepot <- fixed(1.1)
    label("Relative bioavailability of the sprinkle-capsule formulation (theta14, ratio)") # Table 2 theta14 = 1.1 FIXED

    # Covariate effects: Table 2.
    e_wt_vc <- 0.840
    label("Power exponent on (WT / 67.4) for Vc/F (theta8, unitless)") # Table 2 theta8 = 0.840
    e_wt_cl <- 0.669
    label("Power exponent on (WT / 67.4) for CL/F (theta9, unitless)") # Table 2 theta9 = 0.669
    e_age_q <- 0.719
    label("Power exponent on (AGE / 25) for Q/F (theta10, unitless)") # Table 2 theta10 = 0.719
    e_ped9m_ka <- 1.53
    label("Fractional change in ka for patients aged 9 months to <18 years (theta11, unitless)") # Table 2 theta11 = 1.53
    e_dis_all_ka <- -1.67
    label("Exponential effect of the ALL/LL patient type on ka (theta16, unitless)") # Table 2 theta16 = -1.67
    e_dis_all_vc <- -0.564
    label("Exponential effect of the ALL/LL patient type on Vc/F (theta17, unitless)") # Table 2 theta17 = -0.564

    # Non-renal CL/F maturation function taken from the literature (Table 2
    # thetas 12, 13 and 15, all fixed).
    lhill <- fixed(log(0.83))
    label("Shape (Hill) parameter of the non-renal CL/F ontogeny function (theta12, unitless)") # Table 2 theta12 = 0.83 FIXED
    lpna50 <- fixed(log(0.244))
    label("Age at half-maximal non-renal CL/F maturation, AGE50 (theta13, years)") # Table 2 theta13 = 0.244 FIXED
    fbirth_cl <- fixed(0.05)
    label("Fraction of adult non-renal CL/F activity at birth, FBirth (theta15, fraction)") # Table 2 theta15 = 0.05 FIXED

    # Linear PK/PD relationship, Table 4.
    lslope_antiFXa <- log(0.839)
    label("Slope of anti-factor Xa activity vs apixaban plasma concentration (unitless, ng/mL per ng/mL)") # Table 4 Slope = 0.839

    # IIV: Table 2, variances. The ZVC:ZVP row is the Vp/Q covariance
    # (omega 5,4; footnote r = 0.887 for ZQ:ZVP; $OMEGA BLOCK(2) on eta4 V3
    # and eta5 Q in Data S1).
    etalka ~ 0.364 # Table 2 'ZKA' omega 1,1 = 0.364
    etalcl + etalvc ~ c(0.158, 0.101, 0.0766) # Table 2 'ZCL' 0.158, 'ZVC:ZCL' 0.101, 'ZVC' 0.0766
    etalvp + etalq ~ c(1.21, 0.915, 0.880) # Table 2 'ZVP' 1.21, omega 5,4 0.915, 'ZQ' 0.880

    # Residual error: Table 2, additive on log concentration, variance
    # 0.0997 (SD 0.316).
    expSd <- 0.316
    label("Additive residual SD on log-transformed apixaban concentration") # Table 2 sigma 1,1 = 0.0997 (SD 0.316)

    # The PK/PD linear regression residual SD is not reported in the paper.
    addSd_antiFXa <- fixed(0)
    label("Additive residual SD on anti-factor Xa activity, not reported (ng/mL)")
  })
  model({
    # Age-group indicator PED9M: 1 for 9 months to <18 years (Table 2 footnote;
    # control stream 'PED9M=0; IF(PED4.EQ.1)PED9M=1').
    ped9m <- (AGE >= 0.75) * (AGE < 18)

    # Absorption rate constant, Table 2 equations.
    ka <- exp(lka + etalka) * (1 + e_ped9m_ka * ped9m) * exp(e_dis_all_ka * DIS_ALL)

    # Non-renal maturation of CL/F (control stream CLHMAT).
    hill <- exp(lhill)
    pna50 <- exp(lpna50)
    fmat <- fbirth_cl + (1 - fbirth_cl) * AGE^hill / (AGE^hill + pna50^hill)

    cl <- exp(lcl + etalcl) * fmat * (WT / 67.4)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 67.4)^e_wt_vc * exp(e_dis_all_vc * DIS_ALL)
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

    # Dose in mg and volume in L give mg/L; x1000 reports ng/mL (control
    # stream S2 = V2/1000).
    Cc <- 1000 * central / vc

    # Zero-intercept linear PK/PD (Methods, PK/PD analysis: AXA = slope x
    # apixaban concentration); drug-calibrated activity in ng/mL.
    slope_antiFXa <- exp(lslope_antiFXa)
    antiFXa <- slope_antiFXa * Cc

    Cc ~ lnorm(expSd)
    antiFXa ~ add(addSd_antiFXa)
  })
}
