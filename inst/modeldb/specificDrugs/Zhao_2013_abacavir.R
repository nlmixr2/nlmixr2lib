Zhao_2013_abacavir <- function() {
  description <- "Two-compartment population PK model for oral abacavir in HIV-infected infants, toddlers, and children (Zhao 2013); body weight is the only retained covariate (allometric on CL/F and V1/F with estimated exponents and reference weight 17.6 kg)."
  reference <- "Zhao W, Piana C, Danhof M, Burger D, Della Pasqua O, Jacqz-Aigrain E. Population pharmacokinetics of abacavir in infants, toddlers and children. Br J Clin Pharmacol. 2013;75(6):1525-1535. doi:10.1111/bcp.12024"
  vignette <- "Zhao_2013_abacavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling on CL/F and V1/F with estimated exponents (0.802 and 0.810) and reference weight 17.6 kg (population median).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in forward selection on CL/F; produced a significant OFV drop alone but did not survive backward elimination once weight was retained (paper Results section)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a covariate (paper called it 'gender'); not retained in the final model (paper Results)."
    ),
    FORM_TABLET = list(
      description = "Tablet vs solution formulation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested in forward selection on CL/F; no pharmacokinetic differences between the two formulations were retained in the final model (paper Results / Conclusions)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 69,
    n_studies      = 3,
    age_range      = "0.42-12.84 years",
    age_median     = "5.66 years",
    weight_range   = "7.6-60.9 kg",
    weight_median  = "17.6 kg",
    sex_female_pct = NULL,
    race_ethnicity = "Not stratified in source (European and Ugandan paediatric cohorts pooled).",
    disease_state  = "HIV type-1 infected paediatric patients (infants, toddlers, and children).",
    dose_range     = "WHO weight-band paediatric oral abacavir: 16 mg/kg/day weight-normalized or 300 / 450 / 600 mg fixed dose per day (split as 8 mg/kg twice daily or 150-300 mg twice daily; some patients received once-daily dosing). Tablet or oral solution formulation.",
    regions        = "France (PENTA 13), multinational European (PENTA 15), Uganda (ARROW pharmacokinetic substudy).",
    n_observations = 1065,
    notes          = "Pooled meta-analysis of three paediatric clinical studies: PENTA 13 (n = 14, 2-13 years), PENTA 15 (n = 18, 3 months - 3 years), and the ARROW pharmacokinetic substudy (n = 37, 3-12 years). Steady-state pharmacokinetic samples taken at T0, T1, T2, T3, T4, T6, T8, T12 h post-dose (twice daily) plus T24 h for once daily; 138 profiles total. Baseline demographics in Table 1 of the source."
  )

  ini({
    # Structural parameters (Zhao 2013 Table 2 final estimates).
    # Reference weight for allometric scaling is the population median 17.6 kg.
    lka  <- log(0.913); label("Absorption rate constant (ka, 1/h)")                                    # Table 2
    lcl  <- log(20.1);  label("Apparent clearance at 17.6 kg reference (CL/F, L/h)")                   # Table 2
    lvc  <- log(13.0);  label("Apparent central volume of distribution at 17.6 kg reference (V1/F, L)") # Table 2
    lvp  <- log(13.5);  label("Apparent peripheral volume of distribution (V2/F, L)")                  # Table 2
    lq   <- log(2.0);   label("Apparent intercompartmental clearance (Q/F, L/h)")                      # Table 2

    # Allometric exponents (estimated by Zhao 2013, NOT fixed at the canonical 0.75 / 1).
    e_wt_cl <- 0.802; label("Allometric exponent on CL/F (unitless)") # Table 2 (theta_1, RSE 11.6 %)
    e_wt_vc <- 0.810; label("Allometric exponent on V1/F (unitless)") # Table 2 (theta_2, RSE 23.3 %)

    # IIV (omega^2 = log(CV^2 + 1); exponential IIV model per paper Methods).
    # CL/F: BSV (21.9 % CV) and IOV (20.4 % CV) were retained and coupled additively
    # in NONMEM (CL_ij = TVCL * exp(eta_i + kappa_ij)). For single-occasion forward
    # simulation in nlmixr2lib both terms are collapsed onto etalcl by adding the
    # log-scale variances (matches the Birgersson 2016 artemisinin convention).
    # Reproduction: log(0.219^2 + 1) + log(0.204^2 + 1) = 0.04686 + 0.04077 = 0.08763.
    etalcl ~ 0.08763 # Table 2 (combined 21.9 % BSV + 20.4 % IOV on CL/F; see vignette Errata)
    etalvc ~ 0.20509 # Table 2 (47.7 % CV BSV on V1/F)
    etalvp ~ 0.28555 # Table 2 (57.5 % CV BSV on V2/F)
    etalq  ~ 0.16608 # Table 2 (42.5 % CV BSV on Q/F)

    # Residual error (proportional model, paper Methods / Table 2).
    propSd <- 0.382; label("Proportional residual error (fraction)") # Table 2 (38.2 % CV)
  })
  model({
    # Individual PK parameters with allometric weight scaling (reference 17.6 kg).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 17.6)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 17.6)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
