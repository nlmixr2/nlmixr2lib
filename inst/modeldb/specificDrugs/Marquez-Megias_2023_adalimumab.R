`Marquez-Megias_2023_adalimumab` <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and linear elimination for adalimumab in adults with inflammatory bowel disease, with albumin and anti-drug-antibody covariates on apparent clearance (Marquez-Megias 2023)"
  reference <- "Marquez-Megias S, Nalda-Molina R, Mas-Serrano P, Ramon-Lopez A. Population Pharmacokinetic Model of Adalimumab Based on Prior Information Using Real World Data. Biomedicines. 2023;11(10):2822. doi:10.3390/biomedicines11102822"
  vignette <- "Marquez-Megias_2023_adalimumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin (time-varying)",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F using mean-normalization (ALB/3.77)^e_alb_cl per Marquez-Megias 2023 Equation (3) and supplementary Figure S1. Reference 3.77 g/dL is the mean albumin in the studied population (paper text following Equation (3); Table 1 reports the median 3.86 g/dL). Time-varying per Methods (missing values imputed with the patient's own mean).",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-adalimumab antibody (AAA) positivity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AAA-negative)",
      notes              = "Categorical effect on CL/F: CL/F = CL/F_pop * (1 + ADA_POS * e_ada_cl). e_ada_cl = 4.5 was fixed to the reference (Ternant 2015) value because only 9/54 patients (16.7%) were AAA-positive (paper Section 3.2). A patient is AAA-positive if titres exceed 10 ng/mL on at least one occasion (Methods Section 2.2). Source paper labels this covariate 'AAA'; renamed to canonical ADA_POS per inst/references/covariate-columns.md.",
      source_name        = "AAA"
    )
  )

  population <- list(
    n_subjects     = 54L,
    n_studies      = 1L,
    age_range      = "11-89 years",
    age_median     = "43.5 years",
    weight_range   = "34.8-94.0 kg",
    weight_median  = "66.5 kg",
    sex_female_pct = 44.4,
    race_ethnicity = "Not reported (single-center Spanish cohort).",
    disease_state  = "Adults and adolescents with inflammatory bowel disease (Crohn's disease 85.2%, ulcerative colitis 14.8%).",
    dose_range     = "Subcutaneous induction 160/80 mg at weeks 0/2 (43/54 patients) or 80/40 mg (2/54 patients), then maintenance 40 mg every 2 weeks. Induction regimen unknown for 9 patients.",
    regions        = "Single-center retrospective study, Dr. Balmis General University Hospital of Alicante, Spain (2014-2022).",
    bmi_median     = "22.84 kg/m^2 (range 14.1-32.03)",
    albumin_median = "3.86 g/dL (range 1.97-4.96); mean 3.77 g/dL used as the model reference value (mALB).",
    aaa_positive_pct = 16.7,
    originator_pct   = 70.4,
    biosimilar_pct   = 27.8,
    immunomodulator_pct = 40.7,
    notes          = "Therapeutic drug-monitoring dataset comprising 148 trough serum concentrations (19 during induction, 129 during maintenance). Quantification by ELISA (LISA TRACKER Duo Drug + ADAb, TheraDiag); LLOQ 0.1 mg/L for adalimumab and 10 ng/mL for AAA. Final model was developed with informative priors on IIV(CL/F) and IIV(V/F) from the Ternant 2015 reference model (Eur J Clin Pharmacol 71:1155-1157)."
  )

  ini({
    # Structural PK parameters - typical values for the reference patient
    # (mean albumin 3.77 g/dL, AAA-negative). Final model column of
    # Marquez-Megias 2023 Table 3.
    lka <- fixed(log(0.00625)); label("First-order SC absorption rate constant ka (1/h); fixed to Ternant 2015 reference value")  # Marquez-Megias 2023 Table 3 / Section 3.2
    lcl <- log(0.0312);         label("Apparent clearance CL/F for the reference patient (L/h)")                                  # Marquez-Megias 2023 Table 3 (Final model)
    lvc <- log(7.76);           label("Apparent volume of distribution V/F (L)")                                                  # Marquez-Megias 2023 Table 3 (Final model)

    # Covariate effects on CL/F - Equation (3) of Marquez-Megias 2023:
    #   CL/F = CL/F_pop * (1 + AAA * e_ada_cl) * (ALB / 3.77)^e_alb_cl
    # AAA effect (multiplier in 1 + AAA*coef) was fixed to the reference
    # model value because only 9/54 patients were AAA-positive.
    e_ada_cl <- fixed(4.5);     label("AAA multiplicative effect on CL/F (in 1 + AAA*e_ada_cl); fixed to Ternant 2015 reference value") # Marquez-Megias 2023 Section 3.2 (text), Table 3 (Reference Model column reports 4.5)
    e_alb_cl <- -2.33;          label("Power exponent of albumin on CL/F ((ALB/3.77)^e_alb_cl, ALB in g/dL)")                          # Marquez-Megias 2023 Table 3 (Final model)

    # Inter-individual variability. Table 3 reports IIV_CL/F = 0.667 and
    # IIV_V/F = 0.477. Monolix's default population-parameter table prints
    # omega (SD on log scale) for log-normal random effects, so these are
    # interpreted as omega values; the variances passed to nlmixr2 are the
    # squares.
    etalcl ~ 0.667^2  # Marquez-Megias 2023 Table 3 (Final model: omega_CL = 0.667 -> variance 0.4449)
    etalvc ~ 0.477^2  # Marquez-Megias 2023 Table 3 (Final model: omega_V  = 0.477 -> variance 0.2275)

    # Residual error - proportional only. The additive component of the
    # combined error model from the reference model was dropped because its
    # RSE was 83.8% (Section 3.2 of Marquez-Megias 2023).
    propSd <- 0.547;            label("Proportional residual error (fraction)")  # Marquez-Megias 2023 Table 3 (Final model)
  })

  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # Individual PK parameters for the reference patient (mean albumin
    # 3.77 g/dL, AAA-negative). Covariate forms per Equation (3) of
    # Marquez-Megias 2023 / supplementary Figure S1 Monolix code:
    #   CL/F = CL/F_pop * (1 + AAA * e_ada_cl) * (alb_gdL / 3.77)^e_alb_cl
    #   V/F  = V/F_pop
    #   ka   = ka_pop
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + ADA_POS * e_ada_cl) *
      (alb_gdL / 3.77)^e_alb_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, V/F in L -> mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
