Gotta_2014_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 into the central compartment and first-order elimination ",
    "for oral imatinib in European adults with chronic myeloid leukemia ",
    "under routine care, from the largest cohort among the models ",
    "evaluated by Yang 2025 (2478 patients, 4095 samples). CL/F carries a ",
    "sex effect (15.2% lower in women) and a hinged linear age effect that ",
    "peaks at 40 years: clearance rises with age below 40 and falls with ",
    "age above 40. Inter-individual variability is estimated on CL/F and ",
    "Vc/F as a correlated 2x2 block; residual error is combined ",
    "proportional plus additive. TRANSCRIBED FROM A SECONDARY SOURCE: the ",
    "parameter values come from Table 1 of Yang 2025, an external ",
    "evaluation of 15 published imatinib population PK models, not from ",
    "the primary publication. Re-extract from Gotta 2014 when that paper ",
    "is obtained; note that the published inter-individual variability on ",
    "the residual-error magnitude itself is not encodable in nlmixr2 and ",
    "is omitted, as documented in the ini() block."
  )
  reference <- paste0(
    "Gotta V, Bouchet S, Widmer N, Schuld P, Decosterd LA, Buclin T, ",
    "Mahon FX, Csajka C, Molimard M. Large-scale imatinib ",
    "dose-concentration-effect study in CML patients under routine care ",
    "conditions. Leuk Res. 2014;38(7):764-772. ",
    "doi:10.1016/j.leukres.2014.03.023. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Gotta et al. (2014)' ",
    "and Table 1 footnote e. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "Enters CL/F as the multiplicative factor (1 + (-0.152 if ",
        "Female)), i.e. (1 - 0.152 * SEXF) (Yang 2025 Table 1). Women ",
        "therefore have a 15.2% lower typical apparent oral clearance than ",
        "men, and correspondingly higher steady-state trough ",
        "concentrations on the same milligram dose. The source column is ",
        "written 'Sex' in the Covariates-included cell and resolved to the ",
        "female contrast by the equation itself; the canonical SEXF coding ",
        "(1 = female) matches the published contrast directly, so no sign ",
        "flip is needed."
      ),
      source_name        = "Sex"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as the HINGED LINEAR factor (1 + (AGE - 40) * ",
        "theta_Age), where Yang 2025 Table 1 footnote e supplies the two ",
        "branches: 'When Age < 40 the theta_Age = 0.00475, else theta_Age ",
        "= -0.00566 when Age > 40'. Both branches are centred on 40 years, ",
        "so the factor equals exactly 1 at age 40 and the function is ",
        "continuous at the hinge regardless of which branch is taken there ",
        "(the AGE - 40 multiplier is zero); the boundary case is therefore ",
        "immaterial and the model assigns age exactly 40 to the '>= 40' ",
        "branch. Because both slopes multiply a term that is negative ",
        "below the hinge and positive above it, apparent oral clearance ",
        "PEAKS at 40 years and declines in both directions: a 20-year-old ",
        "has a factor of 1 + (-20)(0.00475) = 0.905 and a 70-year-old ",
        "1 + (30)(-0.00566) = 0.830. Note that the factor turns negative ",
        "above age 217, far outside any plausible range, so no guard is ",
        "needed in practice; the model must not, however, be extrapolated ",
        "to children, for whom the cohort provides no support (the study ",
        "age range is 18-91 years)."
      ),
      source_name        = "Age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2478L,
    n_studies      = 1L,
    n_observations = "4095 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "18-91 years",
    disease_state  = "Adults with chronic myeloid leukemia (CML) under routine care conditions",
    dose_range     = "Oral imatinib 100-1200 mg total daily dose",
    regions        = "Europe",
    bioanalytical  = "LC-MS/MS, limit of quantification 25 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "By far the largest cohort among the 15 models evaluated by Yang ",
      "2025 -- 2478 patients against a median of about 60 for the rest -- ",
      "and the widest dose range (100-1200 mg/day). Yang 2025 Results ",
      "describes it as 'founded on observation data collected from 2006 to ",
      "2010'. Demographic detail beyond the row above (weight range, sex ",
      "split, race) is not reported by the secondary source and must be ",
      "read from the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Gotta row) -----
    # The typical value is for the REFERENCE subject: male, aged exactly 40
    # years, at which both covariate factors equal 1.
    lcl <- log(17.3); label("Apparent oral clearance CL/F in a 40-year-old man (L/h)")  # Yang 2025 Table 1: CL/F = 17.3 x (1 + (-0.152 if Female)) x (1 + (Age - 40) x theta_Age)
    lvc <- log(430); label("Central volume of distribution Vc (L)")  # Yang 2025 Table 1: Vc = 430
    ld1 <- log(3.1); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 3.1

    # ----- Covariate effects on CL/F -----
    e_sexf_cl <- -0.152; label("Fractional change in CL/F for female sex (unitless)")  # Yang 2025 Table 1: (1 + (-0.152 if Female))

    # Hinged linear age effect: two slopes about a knot at 40 years.
    e_age_lt40_cl <- 0.00475; label("Slope of (AGE - 40) on CL/F below age 40 (per year)")  # Yang 2025 Table 1 footnote e: When Age < 40 the theta_Age = 0.00475
    e_age_ge40_cl <- -0.00566; label("Slope of (AGE - 40) on CL/F at or above age 40 (per year)")  # Yang 2025 Table 1 footnote e: else theta_Age = -0.00566 when Age > 40

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 37.6%', 'CV%(Vc): 50.4%' and,
    # unusually for this table, the CORRELATION rather than the covariance:
    # 'rho_CL-Vc: 0.73'. The tabulated CV% is taken as omega (the log-scale
    # standard deviation), so the variances are (CV/100)^2 and the
    # covariance is reconstructed as rho * omega_CL * omega_Vc =
    # 0.73 x 0.376 x 0.504 = 0.138338.
    etalcl + etalvc ~ c(0.141376,
                        0.138338, 0.254016)  # Yang 2025 Table 1: omega^2 = 0.376^2 and 0.504^2; Cov = 0.73 x 0.376 x 0.504

    # ----- Residual unexplained variability -----
    # The additive term is tabulated in ng/mL, matching the units of Cc, so
    # no conversion is needed.
    #
    # DOCUMENTED DEVIATION: Yang 2025 Table 1 additionally reports
    # 'CV%(sigma_Prop): 36.2%' for this model -- an inter-individual random
    # effect on the MAGNITUDE of the proportional residual error itself, so
    # that each subject has their own sigma. nlmixr2's error-model syntax
    # takes a fixed residual-SD parameter and provides no supported way to
    # attach an eta to it, so this component is OMITTED here. The
    # consequence is that simulated between-subject spread in residual
    # scatter is narrower than the published model implies; typical-value
    # predictions, the covariate model and the CL/V random effects are all
    # unaffected. This omission is also recorded in the vignette's
    # Assumptions and deviations section.
    propSd <- 0.292; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 29.2%
    addSd <- 82.8; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 82.8 ng/mL
  })

  model({
    # ----- 1. Covariate factors on CL/F -----
    # Hinged linear age effect. age_lt40 selects the below-40 slope; the
    # complement selects the at-or-above-40 slope. The two branches agree
    # at exactly age 40 because the (AGE - 40) multiplier is zero there.
    age_lt40 <- (AGE < 40)
    theta_age <- e_age_lt40_cl * age_lt40 + e_age_ge40_cl * (1 - age_lt40)
    age_factor <- 1 + (AGE - 40) * theta_age

    # Sex effect: SEXF = 1 for women reduces CL/F by 15.2%.
    sex_factor <- 1 + e_sexf_cl * SEXF

    # ----- 2. Individual parameters -----
    cl <- exp(lcl + etalcl) * sex_factor * age_factor
    vc <- exp(lvc + etalvc)
    d1 <- exp(ld1)

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system -----
    # One compartment, zero-order input. Dose records must carry rate = -2
    # so rxode2 uses the modelled duration dur(central) = d1.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # ----- 5. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
