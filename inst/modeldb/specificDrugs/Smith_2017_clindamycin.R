Smith_2017_clindamycin <- function() {
  description <- "One-compartment population PK model for intravenous clindamycin in obese and nonobese children, with allometric total body weight on CL and V, sigmoidal Hill maturation on CL by postmenstrual age, and power effects of serum albumin and alpha-1 acid glycoprotein on V (Smith 2017)."
  reference <- "Smith MJ, Gonzalez D, Goldman JL, Yogev R, Sullivan JE, Reed MD, Anand R, Martz K, Berezny K, Benjamin DK Jr, Smith PB, Cohen-Wolkowiez M, Watt K, on behalf of the Best Pharmaceuticals for Children Act-Pediatric Trials Network Steering Committee. Pharmacokinetics of Clindamycin in Obese and Nonobese Children. Antimicrob Agents Chemother. 2017;61(4):e02014-16. doi:10.1128/AAC.02014-16"
  vignette <- "Smith_2017_clindamycin"
  units <- list(time = "hr", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Reference 70 kg. Allometric scaling: CL ~ (WT/70)^0.75, V ~ (WT/70)^1.0 (both exponents fixed per Smith 2017 Methods Equations 4-5). TBW was found to be the most robust body-size descriptor over FFM/NFM/LBW per Results.",
      source_name        = "TBW"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The Smith 2017 paper reports postmenstrual age in WEEKS in the sigmoidal Hill maturation function on CL (TM50 = 39.5 weeks, Hill = 2.83). Canonical PAGE in nlmixr2lib is in months, so the model converts internally as PMA_weeks = PAGE * 4.345 before evaluating the Hill equation.",
      source_name        = "PMA"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying or baseline (per source dataset). Reference 3.3 g/dL (Smith 2017 Methods Equation 5 typical value used in the V relationship). Power effect on V: (ALB/3.3)^(-0.83).",
      source_name        = "ALB"
    ),
    AAG = list(
      description        = "Alpha-1 acid glycoprotein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying or baseline (per source dataset). Reference 2.4 g/L (= 2.4 mg/mL; Smith 2017 reports AAG in mg/mL, 1 mg/mL = 1 g/L). Power effect on V: (AAG/2.4)^(-0.25).",
      source_name        = "AAG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 220L,
    n_studies      = 3L,
    age_range      = "5 days postnatal (Staph Trio neonates) to 20.5 years (PTN POPS)",
    weight_range   = "0.5-224 kg (pooled across three trials)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in the main paper",
    disease_state  = "Pediatric patients (obese and nonobese) and premature infants receiving intravenous clindamycin per standard of care; indications include skin and soft-tissue infections, bone and joint infections, and other invasive Staphylococcus aureus infections.",
    dose_range     = "30-40 mg/kg/day IV (CLIN01); 10 mg/kg every 6, 8, or 12 h IV (Staph Trio); per standard-of-care dosing (PTN POPS)",
    regions        = "United States (multicenter)",
    notes          = "Pooled cohort of 220 children (76 obese with BMI >= 95th percentile for age) contributing 420 plasma samples across three Best Pharmaceuticals for Children Act-Pediatric Trials Network studies: CLIN01 (n=21 adolescents with BMI >= 85th percentile; NCT01744730), PTN POPS (n=178 pediatric standard-of-care; NCT01431326), and Staph Trio (n=21 neonates; NCT01728363). Demographics per Smith 2017 Table 1; combined cohort summary per Results."
  )

  ini({
    # Structural parameters -- Smith 2017 Table 2 (Final model)
    # Reference subject is a 70 kg adult (CL_70kg, V_70kg) with full maturation (PMA -> infinity)
    # and reference protein binding (ALB = 3.3 g/dL, AAG = 2.4 g/L).
    lcl <- log(13.8); label("Clearance for a 70 kg fully mature adult (CL_70kg, L/hr)")  # Smith 2017 Table 2 (CL_70kg = 13.8 L/h, RSE 6.2%)
    lvc <- log(63.6); label("Central volume of distribution for a 70 kg reference subject (V_70kg, L)")  # Smith 2017 Table 2 (V_70kg = 63.6 L, RSE 5.0%)

    # Allometric exponents on TBW (fixed at canonical 0.75 / 1.0 per Methods Equations 4-5; no RSE reported in Table 2)
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL for TBW (unitless, fixed)")  # Smith 2017 Methods Equation 4 (fixed)
    e_wt_vc <- fixed(1);    label("Allometric exponent on V for TBW (unitless, fixed)")   # Smith 2017 Methods Equation 5 (fixed)

    # Sigmoidal Hill maturation function on CL with postmenstrual age (weeks)
    pma_tm50 <- 39.5; label("Postmenstrual age at 50% mature CL (TM50, weeks)")           # Smith 2017 Table 2 (TM50 = 39.5 weeks, RSE 12.1%)
    pma_hill <- 2.83; label("Hill coefficient for sigmoidal CL maturation (unitless)")    # Smith 2017 Abstract Equation (HILL = 2.83; Table 2 rounded value 2.8)

    # Protein-binding covariate effects on V (power form, reference ALB = 3.3 g/dL, AAG = 2.4 g/L)
    e_alb_vc <- -0.83; label("Power exponent of ALB on V (unitless)")  # Smith 2017 Table 2 (Albumin on V exponent = -0.83, RSE 27.9%)
    e_aag_vc <- -0.25; label("Power exponent of AAG on V (unitless)")  # Smith 2017 Table 2 (Alpha-1 acid glycoprotein on V exponent = -0.25, RSE 44.0%)

    # IIV (correlated CL and V); CV% from Table 2 converted via omega^2 = log(1 + CV^2)
    # IIV(CL) = 58.5% -> log(1 + 0.585^2) = 0.29449
    # IIV(V)  = 11.6% -> log(1 + 0.116^2) = 0.01337
    # Correlation 0.8 -> cov = 0.8 * sqrt(0.29449 * 0.01337) = 0.05020
    etalcl + etalvc ~ c(0.29449,
                        0.05020, 0.01337)  # Smith 2017 Table 2 (IIV CL 58.5%, rho_CL-V 0.8, IIV V 11.6%)

    # Residual error (single representative proportional error; see Assumptions and deviations in the vignette)
    propSd <- 0.336; label("Proportional residual error (fraction)")  # Smith 2017 Table 2 (Prop. PTN POPS = 33.6%; representative value -- see vignette)
  })
  model({
    # Convert canonical PAGE (months) to postmenstrual age in weeks for the sigmoidal Hill
    # maturation equation (Smith 2017 reports TM50 in weeks).
    pma_wk <- PAGE * 4.345

    # Sigmoidal Hill maturation function on CL (Smith 2017 Methods Equation 4)
    maturation_cl <- pma_wk^pma_hill / (pma_tm50^pma_hill + pma_wk^pma_hill)

    # Protein-binding multiplier on V (Smith 2017 Methods Equation 5)
    pb_vc <- (ALB / 3.3)^e_alb_vc * (AAG / 2.4)^e_aag_vc

    # Individual PK parameters with allometric weight scaling (reference 70 kg)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * maturation_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * pb_vc

    kel <- cl / vc

    # IV-only administration (CLIN01, PTN POPS, Staph Trio all dosed IV; no extravascular route).
    d/dt(central) <- -kel * central

    # Concentration: dose in mg, V in L -> mg/L = ug/mL
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
