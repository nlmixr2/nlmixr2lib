Thuo_2011_ciprofloxacin <- function() {
  description <- "One-compartment population PK model with first-order absorption and absorption lag for oral ciprofloxacin in Kenyan children with severe malnutrition (Thuo 2011). Apparent CL and apparent Vc are allometrically scaled to body weight (exponents 0.75 and 1) and modified by linear deviations from a serum sodium reference of 136 mmol/L; apparent CL is further reduced by 28.3% in the paper-defined high-mortality-risk stratum."
  reference   <- "Thuo N, Ungphakorn W, Karisa J, Muchohi S, Muturi A, Kokwaro G, Thomson AH, Maitland K. Dosing regimens of oral ciprofloxacin for children with severe malnutrition: a population pharmacokinetic study with Monte Carlo simulation. J Antimicrob Chemother. 2011 Oct;66(10):2336-45. doi:10.1093/jac/dkr314"
  vignette    <- "Thuo_2011_ciprofloxacin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry in Thuo 2011 (single 48 h dosing window). Used for allometric scaling on apparent CL (exponent 0.75) and apparent Vc (exponent 1) with reference weight 70 kg.",
      source_name        = "WT"
    ),
    SOD = list(
      description        = "Serum sodium concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum sodium at study entry; Thuo 2011 cohort median 136 mmol/L (range 120-160). Enters apparent CL and apparent Vc as linear centered-deviation effects (1 + e * (SOD - 136)).",
      source_name        = "Na+"
    ),
    MORTRISK_HIGH = list(
      description        = "High-mortality-risk composite indicator (Berkley 2003 criteria)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = low or intermediate mortality risk",
      notes              = "Thuo 2011 / Berkley 2003 paediatric severe-malnutrition risk stratum: 1 if depressed conscious state OR bradycardia (HR < 80) OR shock (capillary refill >= 2 s, temperature gradient or weak pulse) OR hypoglycaemia (glucose < 3 mmol/L); 0 otherwise. Intermediate-risk children (deep acidotic breathing, severe dehydration with diarrhoea, lethargy, hyponatraemia, or hypokalaemia) are pooled with low-risk into the 0 reference because the final model only contrasts high-vs-not-high.",
      source_name        = "high risk"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52,
    n_studies      = 1,
    age_range      = "8-102 months",
    age_median     = "23 months (IQR 15-33)",
    weight_range   = "4.1-14.5 kg",
    weight_median  = "6.9 kg (IQR 6.1-8.4)",
    sex_female_pct = 44,
    race_ethnicity = "Kenyan paediatric inpatients (sub-Saharan African); ethnic subdivisions not reported.",
    disease_state  = "Severe malnutrition (weight-for-height z-score <= -3, mid-upper arm circumference < 11 cm, or bilateral pedal oedema); 46% with kwashiorkor; 15% HIV-antibody-positive; risk strata 42% low / 27% intermediate / 31% high mortality risk.",
    dose_range     = "10 mg/kg oral ciprofloxacin every 12 h for 48 h (tablets reformulated into an aqueous suspension by the study pharmacist).",
    regions        = "Kenya (Kilifi District Hospital, Coast Province; KEMRI-Wellcome Trust Research Programme).",
    notes          = "Sparse-sampling popPK study; 202 plasma concentrations from 52 children; baseline median (range): serum sodium 136 (120-160) mmol/L, serum creatinine 44 (27-676) umol/L, estimated CrCl 85.8 (5.0-128.7) mL/min/1.73 m^2 (only one subject with severe renal impairment). Demographics and covariate distributions reproduced from Thuo 2011 Table 1."
  )

  ini({
    # Structural parameters - typical values for a 70 kg adult (allometric reference);
    # Na effects and high-risk effect are centered / multiplicative on the 70 kg typical value.
    lka                <- log(2.97);  label("Absorption rate constant (Ka, 1/h)")                                              # Thuo 2011 Table 2 theta6 = 2.97 /h
    lcl                <- log(42.7);  label("Apparent oral clearance for a 70 kg adult at SOD=136 in low/intermediate-risk children (CL/F, L/h)")  # Thuo 2011 Table 2 theta1 = 42.7 L/h/70 kg
    lvc                <- log(372);   label("Apparent central volume of distribution for a 70 kg adult at SOD=136 (Vc/F, L)")  # Thuo 2011 Table 2 theta4 = 372 L/70 kg
    lalag              <- log(0.742); label("Absorption lag time (Alag, h)")                                                   # Thuo 2011 Table 2 theta7 = 0.742 h

    # Allometric exponents - fixed at standard adult-to-paediatric values per the source paper
    # (Methods 'Pharmacokinetic analysis': CL allometric exponent 0.75; V allometric exponent 1).
    allo_cl <- fixed(0.75); label("Allometric exponent on apparent CL (unitless)")                                            # Thuo 2011 Methods, oral CL allometric equation
    allo_vc <- fixed(1);    label("Allometric exponent on apparent Vc (unitless)")                                            # Thuo 2011 Methods, oral V allometric equation

    # Covariate effects - linear centered-deviation form for sodium; fractional shift for high risk
    e_sod_cl           <-  0.0368; label("Linear sodium effect on apparent CL (per mmol/L from SOD=136)")                     # Thuo 2011 Table 2 theta2 = 0.0368
    e_sod_vc           <-  0.0291; label("Linear sodium effect on apparent Vc (per mmol/L from SOD=136)")                     # Thuo 2011 Table 2 theta5 = 0.0291
    e_mortrisk_high_cl <- -0.283;  label("High-mortality-risk fractional effect on apparent CL")                              # Thuo 2011 Table 2 theta3 = -0.283 (28.3% decrease in CL)

    # IIV (between-subject variability) - log-normal; omega^2 = log(CV^2 + 1)
    # Diagonal block; the source paper investigated covariances but did not retain any.
    etalcl  ~ 0.13554  # 38.1% CV -> log(1 + 0.381^2)
    etalvc  ~ 0.16966  # 43.0% CV -> log(1 + 0.430^2)
    etalka  ~ 0.71304  # 102% CV  -> log(1 + 1.02^2)

    # Residual error - combined additive + proportional in linear concentration space
    propSd <- 0.186;  label("Proportional residual error (fraction)")                                                          # Thuo 2011 Table 2 = 18.6 %CV
    addSd  <- 0.0273; label("Additive residual error (mg/L)")                                                                  # Thuo 2011 Table 2 additive SD
  })

  model({
    # Covariate-effect multipliers
    sod_cl       <- 1 + e_sod_cl * (SOD - 136)
    sod_vc       <- 1 + e_sod_vc * (SOD - 136)
    mortrisk_cl  <- 1 + e_mortrisk_high_cl * MORTRISK_HIGH

    # Individual PK parameters
    ka     <- exp(lka  + etalka)
    cl     <- exp(lcl  + etalcl) * (WT / 70)^allo_cl * sod_cl * mortrisk_cl
    vc     <- exp(lvc  + etalvc) * (WT / 70)^allo_vc * sod_vc
    alag_t <- exp(lalag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag time on the depot compartment
    alag(depot) <- alag_t

    # Concentration: dose in mg, volume in L -> mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
