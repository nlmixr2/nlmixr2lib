Merali_2023_apixaban <- function() {
  description <- "Population PK model for oral apixaban in pooled healthy adults and paediatric subjects from birth to <18 years at risk of venous or arterial thrombosis (Merali 2023): two-compartment with first-order absorption and elimination, allometric body weight on CL/F and Vc/F, a fixed CYP3A4 maturation function of age on CL/F, age on Q/F, a 9-month-to-18-year shift in ka, dose-dependent relative bioavailability above 2.5 mg and a neonatal-formulation bioavailability factor, plus a linear anti-factor Xa activity relationship to apixaban concentration."
  reference <- "Merali SJ, Byon W, Patel YT, Elsrougy A, Marchisin D, Perera V, Chen W, He B, Murthy B. Evaluation of safety, pharmacokinetics, and pharmacodynamics of apixaban in pediatric subjects at risk of venous or arterial thrombotic disorder. CPT Pharmacometrics Syst Pharmacol. 2023;12(4):500-512. doi:10.1002/psp4.12935"
  vignette <- "Merali_2023_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (exponent 0.71) and Vc/F (exponent 0.879) with reference weight 67.4 kg (paper equations for CL_i/F_i and VC_i/F_i). Time-fixed at baseline ('BBWT' in the paper).",
      source_name = "BBWT"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives three effects: (1) the fixed CYP3A4 maturation function on CL/F (age at half-maturation 0.244 years); (2) a power effect on Q/F with reference age 25 years; (3) the 9-month-to-<18-year indicator PED9M on ka, derived in model() as (AGE >= 0.75) * (AGE < 18). Postnatal age in years; the neonate enrolled was 15 days old at dosing.",
      source_name = "AGE"
    ),
    DOSE_APIXABAN_MG = list(
      description = "Administered apixaban dose amount",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent reduction in relative bioavailability above 2.5 mg (TVF1 equations). F1 is 1 at doses of 2.5 mg or less; the 47.5 mg denominator scales the excess to 1 at the 50 mg maximum administered dose. Paediatric single doses were 0.1 to about 3.2 mg; the pooled adult studies dosed 2.5 to 50 mg.",
      source_name = "DOSE"
    ),
    FORM_CAPSULE = list(
      description = "Neonatal apixaban 0.1 mg sprinkle capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (apixaban 0.4 mg/mL oral solution in subjects 28 days to <18 years, or adult tablet)",
      notes = "Multiplies relative bioavailability by the fixed factor theta14 = 1.1 ('F1 for neonatal formulation, which is equal to 1.1 for subjects receiving neonatal formulation and 1 for others'). The neonatal formulation was a 0.1 mg sprinkle capsule dispersed in breast milk or formula, given only to the birth-to-27-day group.",
      source_name = "neonatal formulation"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 254,
    n_studies = 10,
    n_observations = 5220,
    age_range = "15 days to 18 years (paediatric); adult study means 26.7 to 42.0 years",
    age_mean = "Paediatric study CV185118: 5.4 years (SD 5.4); CV185079: 13.9 years (SD 3.9)",
    weight_range = "Paediatric group means 3.1 kg (neonate) to 49.9 kg (12 to <18 years); adult study means 63.1 to 95.6 kg",
    sex_female_pct = 57,
    race_ethnicity = c(White = 86, Black = 10, Other = 4),
    disease_state = "Healthy adults (and two renally impaired adults) pooled with paediatric subjects at risk of venous or arterial thrombotic disorder (CV185118) or with an indwelling central venous catheter (CV185079)",
    dose_range = "Paediatric single oral dose 1.08 to 2.43 mg/m^2 (0.1 mg sprinkle capsule in neonates; 0.4 mg/mL solution otherwise); adults 2.5 to 50 mg single or 7-day multiple doses",
    regions = "United States and international sites",
    notes = "Pooled PK dataset of eight adult phase I studies (198 subjects) and two paediatric phase I studies (56 subjects), Table S2. The sex and race percentages are for the 49 dosed paediatric subjects of CV185118 (Table 1). Reference subject for the covariate effects: 67.4 kg, age 25 years for Q/F, fully mature CYP3A4, single dose of 2.5 mg or less of a non-neonatal formulation."
  )

  ini({
    # Structural parameters: Table 2.
    lka <- log(0.502)
    label("Absorption rate constant for adults and children under 9 months (ka, 1/h)") # Table 2 theta1 = 0.502
    lcl <- log(4.7)
    label("Apparent clearance for a 67.4-kg fully mature subject (CL/F, L/h)") # Table 2 theta2 = 4.7
    lvc <- log(31.3)
    label("Apparent central volume of distribution for a 67.4-kg subject (Vc/F, L)") # Table 2 theta3 = 31.3
    lq <- log(1.4)
    label("Apparent intercompartmental clearance at age 25 years (Q/F, L/h)") # Table 2 theta4 = 1.4
    lvp <- log(18)
    label("Apparent peripheral volume of distribution (Vp/F, L)") # Table 2 theta5 = 18

    # Covariate effects: Table 2.
    e_wt_vc <- 0.879
    label("Power exponent on (WT / 67.4) for Vc/F (unitless)") # Table 2 theta8 = 0.879
    e_wt_cl <- 0.71
    label("Power exponent on (WT / 67.4) for CL/F (unitless)") # Table 2 theta9 = 0.71
    e_age_q <- 0.663
    label("Power exponent on (AGE / 25) for Q/F (unitless)") # Table 2 theta10 = 0.663
    e_ped9m_ka <- 1.42
    label("Proportional shift in ka for subjects aged 9 months to <18 years (unitless)") # Table 2 theta11 = 1.42

    # Fixed CYP3A4 maturation function on CL/F (Bouillon-Pichault 2011,
    # reference 10 of the paper): fraction of adult CL/F =
    # birth + (1 - birth) * AGE^hill / (AGE^hill + t50^hill).
    cyp3a4_hill <- fixed(0.83)
    label("CYP3A4 maturation sigmoidicity exponent (unitless)") # Table 2 theta12 = 0.83 FIXED
    cyp3a4_t50 <- fixed(0.244)
    label("Age at half-maximal CYP3A4 maturation (years)") # Table 2 theta13 = 0.244 FIXED
    cyp3a4_birth <- fixed(0.05)
    label("Fraction of adult CL/F at birth (unitless)") # Table 2 theta15 = 0.05 FIXED

    # Relative bioavailability: dose-dependent reduction above 2.5 mg and a
    # neonatal-formulation factor (paper equations for TVF1 and F1).
    lgamma <- log(0.884)
    label("Shape parameter of the dose-relative-bioavailability relationship (unitless)") # Table 2 theta6 = 0.884
    logitifdepot50 <- -0.465
    label("Logit of the maximum fractional reduction in relative bioavailability, reached at a 50-mg dose (unitless)") # Table 2 theta7 = -0.465
    fdepot_capsule <- fixed(1.1)
    label("Relative bioavailability of the neonatal sprinkle-capsule formulation (unitless)") # Table 2 theta14 = 1.1 FIXED

    # Anti-factor Xa activity: linear regression on apixaban concentration
    # (Discussion; Figure 4). The slope is the paper's overall paediatric
    # estimate; the intercept is not printed and was digitised from the
    # regression line drawn in Figure 4.
    lrbase_antiFXa <- log(0.074)
    label("Intercept of the anti-factor Xa activity regression (IU/mL)") # Digitised from Figure 4 regression line (fitted slope 0.0155 reproduces the printed slope)
    lslope_antiFXa <- log(0.0155)
    label("Slope of anti-factor Xa activity on apixaban concentration (IU/mL per ng/mL)") # Discussion 'overall slope estimate of 0.0155 IU/ng'

    # IIV: Table 2, variances. ka is independent; CL, Vc, Vp and Q share a
    # full covariance block (lower triangle by row: CL; Vc:CL, Vc; Vp:CL,
    # Vp:Vc, Vp; Q:CL, Q:Vc, Q:Vp, Q). The printed covariances reproduce the
    # Table 2 footnote r^2 values, but at their printed precision the block
    # is slightly indefinite (smallest eigenvalue -2.4e-5), so the
    # off-diagonals are scaled by 0.99; the variances are as printed.
    etalka ~ 0.296
    etalcl + etalvc + etalvp + etalq ~ c(
      0.0864,
      0.99 * 0.0468, 0.034,
      0.99 * 0.0789, 0.99 * 0.05, 0.525,
      0.99 * 0.105, 0.99 * 0.0924, 0.99 * 0.393, 0.433
    )

    # Residual error: Table 2 'Exponential' sigma = 0.0844 (variance), SD 0.29.
    expSd <- 0.2905
    label("Exponential (log-scale) residual error SD (unitless)") # Table 2 sigma1,1 = 0.0844 variance, sqrt = 0.2905
  })
  model({
    # Indicator for ages 9 months to <18 years (PED9M); 0 below 9 months and
    # in adults.
    ped9m <- (AGE >= 0.75) * (AGE < 18)
    ka <- exp(lka + etalka) * (1 + ped9m * e_ped9m_ka)

    # Fixed CYP3A4 maturation of CL/F.
    fmat <- cyp3a4_birth +
      (1 - cyp3a4_birth) * AGE^cyp3a4_hill / (AGE^cyp3a4_hill + cyp3a4_t50^cyp3a4_hill)
    cl <- exp(lcl + etalcl) * (WT / 67.4)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 67.4)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq) * (AGE / 25)^e_age_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Relative bioavailability: 1 at doses of 2.5 mg or less; above 2.5 mg
    # reduced by up to expit(logitifdepot50) at 50 mg. The excess is clamped
    # at zero so lower doses cannot raise a negative base to a power.
    gamma <- exp(lgamma)
    ifdepot50 <- expit(logitifdepot50)
    doseExcess <- (DOSE_APIXABAN_MG > 2.5) * (DOSE_APIXABAN_MG - 2.5)
    frel <- 1 - ifdepot50 * (doseExcess / 47.5)^gamma
    fform <- 1 + FORM_CAPSULE * (fdepot_capsule - 1)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel * fform

    # Amounts in mg, volumes in L: central / vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc

    # Anti-factor Xa activity (heparin-calibrated assay, IU/mL); typical-value
    # relationship only, as the paper reports no variability for it.
    rbase_antiFXa <- exp(lrbase_antiFXa)
    slope_antiFXa <- exp(lslope_antiFXa)
    antiFXa <- rbase_antiFXa + slope_antiFXa * Cc

    Cc ~ lnorm(expSd)
  })
}
