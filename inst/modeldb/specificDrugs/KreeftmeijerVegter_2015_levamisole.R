KreeftmeijerVegter_2015_levamisole <- function() {
  description <- paste(
    "One-compartment oral PK model for levamisole in 38 children with",
    "steroid-sensitive nephrotic syndrome (Kreeftmeijer-Vegter 2015,",
    "EudraCT 2005-005745-18). First-order absorption, first-order",
    "elimination, allometric scaling of CL/F (exponent 0.75) and V/F",
    "(exponent 1) to 70 kg, and a linear proportional age effect on",
    "CL/F centred on the population median age of 6.28 years",
    "(-10.1% per additional year). The typical ka (1.2 1/h) was fixed",
    "in the final model with IIV retained. IIV on V/F was modelled as",
    "perfectly correlated with IIV on CL/F (single eta scaled to V/F),",
    "encoded here as a full omega block with covariance equal to",
    "sqrt(var_CL * var_V).",
    sep = " "
  )
  reference <- paste(
    "Kreeftmeijer-Vegter AR, Dorlo TPC, Gruppen MP, de Boer A,",
    "de Vries PJ (2015).",
    "Population pharmacokinetics of levamisole in children with",
    "steroid-sensitive nephrotic syndrome.",
    "British Journal of Clinical Pharmacology 79(6):970-977.",
    "doi:10.1111/bcp.12607.",
    sep = " "
  )
  vignette <- "KreeftmeijerVegter_2015_levamisole"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor; reference weight 70 kg per",
        "Methods (PK data analysis paragraph) and Table 3 row 2",
        "(allometric scaling of CL/F at exponent 0.75 and V/F at",
        "exponent 1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (chronological time since birth)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear proportional effect on CL/F centred on the",
        "population median age of 6.28 years (Table 2 / final-model",
        "equation in Results, Section 'A stepwise covariate analysis').",
        "Slope -0.101 per year, i.e. -10.1% change in CL/F per",
        "additional life year above the median."
      ),
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 1L,
    age_range      = "2.35-13.10 years",
    age_median     = "6.28 years (Table 1)",
    weight_range   = "11-68 kg",
    weight_median  = "21 kg (Table 1)",
    sex_female_pct = 100 * 11 / 38,
    race_ethnicity = c(Caucasian = 47.4, Asian = 50.0, Unknown = 2.6),
    disease_state  = paste(
      "Children with frequently relapsing idiopathic",
      "steroid-sensitive nephrotic syndrome (SSNS), either with or",
      "without steroid dependency."
    ),
    dose_range     = paste(
      "Oral levamisole 2.5 mg/kg every other day (maximum 150 mg),",
      "delivered as 5, 10, 25 or 50 mg film-coated tablets per the",
      "weight-band dosing schedule of Kreeftmeijer-Vegter 2014;",
      "mean administered dose 2.45 mg/kg (SD 0.24).",
      "Treatment continued for 12 months (cohort phase up to a",
      "further 12 months for responders)."
    ),
    regions        = paste(
      "India 47.4%, Netherlands 15.8%, Belgium 13.2%,",
      "France 10.5%, Poland 10.5%, Italy 2.6% (Table 1)."
    ),
    n_observations = 136L,
    sampling       = paste(
      "Sparse sampling at four visits (weeks 8, 12, 20, 24);",
      "in weeks 8 and 20 one pre-dose and one post-dose sample at a",
      "computer-randomised 1, 2, 4 or 6 h post-dose time; in weeks 12",
      "and 24 one sample on a medication-free day (24 h post-dose).",
      "Pre-dose and 24 h post-dose samples (102 of 121) fell below",
      "the assay LOQ and were excluded; final dataset 136 samples."
    ),
    notes          = paste(
      "Steady-state every-other-day dosing in routine clinical care;",
      "sampling restricted to the absorption / peak phase",
      "(0-6 h post-dose) because of paediatric blood-volume",
      "constraints. Storage-time correction (0.0005 per day at",
      "-20 C) was applied directly to the observed concentrations",
      "before NONMEM fitting and is NOT encoded as a model term",
      "(see vignette Errata)."
    )
  )

  ini({
    # Structural parameters -- reference values for 70 kg body weight and
    # the population median age of 6.28 years. From Table 2 of the
    # source.
    lka <- fixed(log(1.2))
    label("Absorption rate constant ka (1/h; fixed in final model)")  # Table 2: ka = 1.2 1/h fixed; pre-covariate point estimate carried into the final model per Results paragraph 3
    lcl <- log(44)
    label("Apparent clearance CL/F at 70 kg, age 6.28 y (L/h)")        # Table 2: CL/F = 44 L/h/70 kg (RSE 8.5%)
    lvc <- log(236)
    label("Apparent central volume V/F at 70 kg (L)")                  # Table 2: V/F = 236 L/70 kg (RSE 13.3%)

    # Allometric exponents -- paper Methods / Results "Allometric scaling
    # ... (power value of 0.75) to standard body weight (70 kg) ... than
    # linear scaling (corresponding allometric power of 1) of CL/F and
    # V/F." Both held fixed at the standard Anderson-Holford values.
    allo_cl <- fixed(0.75)
    label("Allometric exponent on CL/F (unitless; fixed)")             # Results paragraph "Allometric scaling of CL/F (power value of 0.75)"
    allo_vc <- fixed(1)
    label("Allometric exponent on V/F (unitless; fixed)")              # Table 3 row 2 "Allometric scaling of CL/F and V/F to BW (70 kg)"; linear scaling on V/F

    # Linear proportional age effect on CL/F, centred at the population
    # median age (6.28 y). Final-model equation in Results:
    #   CL/F_i = theta_CL/F * (1 + e_age_cl * (AGE - 6.28)) * (WT/70)^0.75
    # The published equation reproduces as
    #   theta_CL/F * -0.101 * (AGE - 6.28) * (WT/70)^0.75
    # but this is a typographical compression of the additive form
    # ("proportional decrease in CL/F per life year"; Results
    # paragraph 4, "the model shows a proportional decrease in CL/F
    # with increasing age in addition to the allometrically scaled
    # effect of body weight").
    e_age_cl <- -0.101
    label("Linear-deviation slope of AGE on CL/F (per year)")          # Table 2: -10.1% per life year (RSE 25%)
    age_ref  <- fixed(6.28)
    label("Reference age for CL/F centring (years)")                   # Table 2 caption: population median age 6.28 y

    # Inter-individual variability. The paper reports IIV as %CV; the
    # internal omega^2 = log(1 + CV^2) under a log-normal IIV model.
    # Table 2 values:
    #   CL/F: CV 31.6% -> omega^2 = log(1 + 0.316^2) = 0.09531
    #   V/F : CV 41.7% -> omega^2 = log(1 + 0.417^2) = 0.16030
    #   ka  : CV 92.2% -> omega^2 = log(1 + 0.922^2) = 0.61509
    # Methods (interindividual variability paragraph): "Because the
    # interindividual variability of V/F correlated 100% with the
    # interindividual variability of CL/F (resulting in
    # overparameterization of our model), the interindividual
    # variability for CL/F was estimated and used to estimate V/F
    # using an additional scaling parameter." Encoded here as a full
    # omega block with rho fixed at 0.99 (rather than exactly 1.0)
    # so the omega matrix is positive definite and rxode2's Cholesky
    # decomposition succeeds during stochastic simulation; the paper's
    # singular rho = 1 form is structurally equivalent to a single
    # eta on CL scaled to V. See vignette Errata.
    # cov_pd = 0.99 * sqrt(0.0953 * 0.1603) = 0.99 * 0.12354 = 0.1223
    etalcl + etalvc ~ c(0.0953,
                        0.1223, 0.1603)                                # Table 2 (rho clipped to 0.99 from 1.0 for PD omega -- see vignette Errata)
    etalka          ~ 0.6151                                           # Table 2: ka CV 92.2%

    # Residual error -- proportional error model (Methods, "Residual
    # variability was modelled using a proportional error model").
    propSd <- 0.207
    label("Proportional residual error (fraction)")                    # Table 2: 20.7% (RSE 25.4%)
  })

  model({
    # Linear proportional age effect on CL/F, centred at age_ref.
    age_eff_cl <- 1 + e_age_cl * (AGE - age_ref)

    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl * age_eff_cl
    # V/F shares the CL eta perfectly (paper Methods); the block omega
    # above carries the correlation, so etalvc is used directly here.
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_vc

    kel <- cl / vc

    # ODE system: 1-compartment with first-order absorption from depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability: F is not separately identifiable (oral dose,
    # parameters reported as CL/F and V/F); rxode2 default F = 1 on the
    # depot is consistent with the published parameterisation.

    # Dose in mg, V/F in L => central/vc in mg/L; convert to ng/mL by
    # multiplying by 1000.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
