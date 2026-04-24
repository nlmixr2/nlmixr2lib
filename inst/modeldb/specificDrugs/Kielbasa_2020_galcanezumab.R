Kielbasa_2020_galcanezumab <- function() {
  description <- "One-compartment population PK model for galcanezumab (humanized IgG anti-CGRP mAb) with first-order SC absorption, linear elimination, and allometric body weight scaling on CL/F (Kielbasa 2020)"
  reference <- "Kielbasa W, Quinlan T. Population Pharmacokinetics of Galcanezumab, an Anti-CGRP Antibody, Following Subcutaneous Dosing to Healthy Individuals and Patients With Migraine. J Clin Pharmacol. 2020;60(2):229-239. doi:10.1002/jcph.1511"
  vignette <- "Kielbasa_2020_galcanezumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric (power) scaling on CL/F with exponent 0.601 and reference weight 73.6 kg (population median). Body weight had a non-qualifying effect on V/F in the source analysis and is not used to scale V/F in the final model.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects        = 1889,
    n_studies         = 7,
    age_range         = "17-65 years",
    age_mean          = "41 years (SD 12)",
    weight_range      = "40-135.5 kg",
    weight_mean       = "75.8 kg (SD 16.8)",
    weight_median     = "73.6 kg",
    sex_female_pct    = 80,
    race_ethnicity    = c(White = 72, NonHispanic = 74),
    disease_state     = "Healthy adults and adults with episodic or chronic migraine (60% episodic, 29% chronic, 11% healthy).",
    dose_range        = "5-300 mg SC, single dose or Q4W / QM regimens; phase 3 studies used 240 mg loading dose then 120 mg QM or 240 mg QM.",
    regions           = "Pooled analysis across 7 clinical studies (phase 1-3); regional breakdown not reported.",
    analysis_dataset  = "15770 PK observations from 1889 individuals included in model development; 1448 observations from 270 patients in Study I5Q-MC-CGAJ used as external validation.",
    notes             = "Baseline demographics from Kielbasa 2020 Results paragraph 1 and Supplemental Table S1."
  )

  ini({
    # Structural parameters - reference values at population median body weight 73.6 kg.
    # Source values are reported on an hourly time base (ka in 1/h, CL/F in L/h);
    # converted to a daily time base for consistency with other mAb models in nlmixr2lib.
    # ka:    0.0199 /h * 24 h/day = 0.4776 /day
    # CL/F:  0.00785 L/h * 24 h/day = 0.1884 L/day
    # V/F:   7.33 L (unchanged)
    lka <- log(0.4776); label("Absorption rate constant (ka, 1/day)")                       # Kielbasa 2020 Table 3 (0.0199 /h * 24)
    lcl <- log(0.1884); label("Apparent clearance for a 73.6 kg subject (CL/F, L/day)")     # Kielbasa 2020 Table 3 (0.00785 L/h * 24)
    lvc <- log(7.33);   label("Apparent volume of distribution (V/F, L)")                   # Kielbasa 2020 Table 3

    # Allometric (power-model) exponent for body weight on CL/F.
    # Model form (equation 2 of Kielbasa 2020): CL/F = theta1 * (WT/MED)^theta2
    # with MED = 73.6 kg (population median). Body weight was tested on V/F as well
    # but did not decrease V/F IIV by > 5% and was not retained in the final model.
    allo_cl_wt <- 0.601; label("Allometric exponent on CL/F for body weight (unitless)")    # Kielbasa 2020 Table 3

    # Inter-individual variability - full 3x3 omega block on ka, CL/F, V/F.
    # The source paper reports two covariances (ka/CL/F and CL/F/V/F); the ka-V/F
    # covariance was not estimated and is fixed to zero. Resulting matrix is PD.
    # omega^2 = log(CV^2 + 1) for log-normally distributed IIV:
    #   ka  92% CV -> log(0.92^2 + 1) = 0.61309
    #   CL  34% CV -> log(0.34^2 + 1) = 0.10942
    #   V   34% CV -> log(0.34^2 + 1) = 0.10942
    # Reported covariances on the omega scale: 0.0694 (ka,CL/F); 0.0716 (CL/F,V/F).
    etalka + etalcl + etalvc ~ c(0.61309,
                                 0.0694,  0.10942,
                                 0,       0.0716,  0.10942)                                 # Kielbasa 2020 Table 3

    # Residual error - proportional only.
    propSd <- 0.22; label("Proportional residual error (fraction)")                         # Kielbasa 2020 Table 3
  })
  model({
    # Individual PK parameters with allometric body-weight scaling on CL/F only.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 73.6)^allo_cl_wt
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
