Archary_2019_lamivudine <- function() {
  description <- "One-compartment population PK model for lamivudine in severely malnourished HIV-infected children (Archary 2019); CL/F matures with age via a sigmoid Emax function, Vc/F decreases linearly with serum triglyceride, and ka steps up between day 1 and day 14 of antiretroviral treatment"
  reference <- "Archary M, McIlleron H, Bobat R, LaRussa P, Sibaya T, Wiesner L, Hennig S. Population pharmacokinetics of abacavir and lamivudine in severely malnourished human immunodeficiency virus-infected children in relation to treatment outcomes. Br J Clin Pharmacol. 2019;85(8):1881-1890. doi:10.1111/bcp.13998"
  vignette <- "Archary_2019_lamivudine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling on CL/F (exponent 0.75) and Vc/F (exponent 1) with reference weight 7 kg (population median).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (chronological time since birth)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; drives the CL/F maturation function (sigmoid Emax with Hill coefficient 1.47 and TM50 = 0.25 y).",
      source_name        = "AGE"
    ),
    TRIG = list(
      description        = "Serum triglyceride concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying or baseline (paper does not separate). Linear-deviation effect on Vc/F centred at 5.3 mmol/L (paper-reported cohort average per the equation; cohort median in Table 1 is 2.2-2.3 mmol/L -- see Errata).",
      source_name        = "TRIG"
    ),
    DAY14 = list(
      description        = "Day-14-post-ART landmark indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (day 1 / pre-rehabilitation)",
      notes              = "Within-subject step indicator: 0 on day 1 of antiretroviral treatment (acute malnutrition baseline), 1 on day 14 (post-nutritional-rehabilitation steady state). Gates the day-1 (0.30) vs day-14 (0.34) typical ka step reported in Table 2 of the source.",
      source_name        = "DAY14"
    )
  )

  population <- list(
    n_subjects     = 75,
    n_studies      = 1,
    age_range      = "0.1-10.8 years (median 1.4)",
    age_median     = "1.4 years",
    weight_range   = "1.88-19.6 kg",
    weight_median  = "7 kg",
    sex_female_pct = 41,
    race_ethnicity = "South African paediatric cohort (race not stratified in source).",
    disease_state  = "Severely malnourished HIV-infected children (weight-for-length Z-scores < -3, mid-upper arm circumference < 115 mm, or peripheral oedema) initiating antiretroviral treatment.",
    dose_range     = "WHO weight-band paediatric oral lamivudine + abacavir + LPV/r dosing; liquid formulation for children < 14 kg, solid formulation > 14 kg (only 2 patients received solid formulations).",
    regions        = "South Africa (King Edward VIII Hospital, Durban).",
    n_observations = 627,
    notes          = "MATCH (Malnutrition and ART Timing in Children with HIV) trial, PACTR201609001751384; 75 patients with day-1 abacavir + lamivudine concentrations, 69 of whom had day-14 samples; 627 lamivudine concentrations sampled 0.8-12.4 h post-dose. Demographics summarised in Table 1 of the source."
  )

  ini({
    # Structural parameters; reference weight 7 kg = population median
    lka  <- log(0.30);   label("Absorption rate constant on day 1 (ka, 1/h)")                   # Table 2 (day-1 typical value)
    lcl  <- log(12.2);   label("Apparent clearance at 7 kg reference, fully matured (CL/F, L/h)") # Table 2
    lvc  <- log(8.22);   label("Apparent central volume at 7 kg reference, TRIG = 5.3 mmol/L (Vc/F, L)") # Table 2

    # Allometric exponents (paper-fixed per Methods Section 2.3)
    e_wt_cl <- 0.75; label("Allometric exponent on CL/F (unitless; fixed)")                     # Methods 2.3
    e_wt_vc <- 1;    label("Allometric exponent on Vc/F (unitless; fixed)")                     # Methods 2.3

    # Maturation function for CL/F (sigmoid Emax in age)
    mat_hill   <- 1.47; label("Hill coefficient of the CL/F sigmoid Emax age-maturation function (unitless)") # Table 2
    mat_age50  <- 0.25; label("Age at 50% mature CL/F (years; ~4 months)")                                   # Table 2

    # Triglyceride effect on Vc/F (linear deviation, ref 5.3 mmol/L)
    e_trig_vc  <- -0.13; label("Linear-deviation slope of TRIG on Vc/F (per mmol/L)")            # Section 3.3 / Table 2 (-13.3% per +1 mmol/L)
    trig_ref   <- 5.3;   label("Reference triglyceride concentration for Vc/F centring (mmol/L)") # Section 3.3

    # DAY14 effect on ka (multiplicative additive shift; encodes the day-1 -> day-14 step)
    # 0.34 / 0.30 - 1 = 0.1333
    e_day14_ka <- 0.1333; label("Effect of DAY14 (post-ART rehabilitation) on ka (fraction)")    # Table 2 (day-14 0.34 vs day-1 0.30 /h)

    # IIV (omega^2 = log(CV^2 + 1))
    # Block omega for CL/F and Vc/F with correlation rho = 0.674
    # var(CL) = log(1 + 0.423^2) = 0.16459
    # var(Vc) = log(1 + 0.66^2)  = 0.36157
    # cov     = rho * sqrt(var(CL) * var(Vc)) = 0.674 * sqrt(0.16459 * 0.36157) = 0.16456
    etalcl + etalvc ~ c(0.1646,
                        0.1646, 0.3616)  # Table 2; rho = 0.674 reported as 'Covariance IIV CL/F and IIV Vc/F'
    # Source reports IOV CL/F = 59.7%; the model file does not encode IOV (would require
    # per-occasion eta multiplexing, see Jonsson_2011_ethambutol.R) -- documented in vignette Errata.

    # Residual error
    propSd <- 0.36; label("Proportional residual error (fraction)")                              # Table 2
  })
  model({
    # Maturation function for CL/F (sigmoid Emax in age, approaches 1 at full maturity)
    maturation_cl <- AGE^mat_hill / (mat_age50^mat_hill + AGE^mat_hill)

    # Triglyceride effect on Vc/F (linear deviation around trig_ref)
    trig_vc <- 1 + e_trig_vc * (TRIG - trig_ref)

    # PK parameters (allometric weight scaling, day-14 ka step)
    ka <- exp(lka) * (1 + e_day14_ka * DAY14)
    cl <- exp(lcl + etalcl) * (WT / 7)^e_wt_cl * maturation_cl
    vc <- exp(lvc + etalvc) * (WT / 7)^e_wt_vc * trig_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, volume in L -> mg/L (= ug/mL)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
