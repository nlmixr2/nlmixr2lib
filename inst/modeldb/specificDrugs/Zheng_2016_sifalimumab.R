Zheng_2016_sifalimumab <- function() {
  description <- "Two-compartment population PK model for sifalimumab (anti-IFN-alpha human IgG1 monoclonal antibody) in adults with systemic lupus erythematosus following repeat fixed intravenous doses (Zheng 2016)."
  reference <- "Zheng B, Yu X-Q, Greth W, Robbie GJ. Population pharmacokinetic analysis of sifalimumab from a clinical phase IIb trial in systemic lupus erythematosus patients. Br J Clin Pharmacol. 2016;81(5):918-928. doi:10.1111/bcp.12864"
  vignette <- "Zheng_2016_sifalimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for power-form covariate scaling of CL and V1 with reference weight 64.3 kg (median of the phase IIb SLE cohort). Exponents estimated (not fixed) at 0.45 on CL and 0.36 on V1 (Zheng 2016 Table 2).",
      source_name        = "WT"
    ),
    BGENE21 = list(
      description        = "Baseline 21-gene type I interferon signature score (whole-blood transcriptomic composite)",
      units              = "fold-change (relative to healthy-donor reference)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for power-form scaling of CL with reference value 12.04 (median of the SLE cohort; range 0.32-38.59). The 21-gene panel is the AstraZeneca/MedImmune SLE development programme signature; effect exponent 0.09 on CL (Zheng 2016 Table 2).",
      source_name        = "BGENE21"
    ),
    STEROID_BL = list(
      description        = "Baseline (concomitant) systemic corticosteroid use at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on steroids at baseline)",
      notes              = "85% of the phase IIb cohort were on baseline steroids. Applied as a multiplicative fractional effect on CL `(1 + 0.11 * STEROID_BL)` and on V1 `(1 - 0.09 * STEROID_BL)` (Zheng 2016 Table 2).",
      source_name        = "BSTEROID"
    ),
    DOSE = list(
      description        = "Subject's assigned fixed dose level (200, 600, or 1200 mg q4w IV)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Treated as a per-subject covariate that takes the same value across all observations for a given subject. Power-form scaling of V1 with reference dose 600 mg and exponent 0.06 (Zheng 2016 Table 2). The effect is small and included here to reproduce the published final model.",
      source_name        = "Dose"
    )
  )

  population <- list(
    n_subjects         = 298L,
    n_studies          = 1L,
    age_range          = "18-73 years",
    age_median         = "40 years",
    weight_range       = "39-131.3 kg",
    weight_median      = "64.3 kg",
    sex_female_pct     = 92,
    race_ethnicity     = c(White = 59, Asian = 15, `African American` = 7, `Native American` = 4, Other = 15),
    disease_state      = "Adults with moderate-to-severe systemic lupus erythematosus (SLE) meeting ACR criteria.",
    dose_range         = "200, 600, or 1200 mg as a 30-60 minute IV infusion every 4 weeks for 52 weeks (13 q4w doses) with an additional loading dose on day 15; 14 IV doses total.",
    regions            = "United States 69%; other countries 31% (multi-region phase IIb trial).",
    baseline_steroid_use_pct = 85,
    bgene21_range      = "0.32-38.59 (median 12.04)",
    sledai_range       = "6-27 (median 10)",
    n_observations     = 3961L,
    notes              = "Baseline demographics per Zheng 2016 Table 1. PK dataset: 3961 quantifiable serum concentrations from 298 sifalimumab-treated patients (102 at 200 mg, 102 at 600 mg, 94 at 1200 mg) out of 431 enrolled in the phase IIb SLE trial (MI-CP180)."
  )

  ini({
    # Structural parameters - Zheng 2016 Table 2 final-model estimates for a typical patient
    # representing median WT = 64.3 kg, median BGENE21 = 12.04, 600 mg dose (for V1), and
    # not on baseline steroids (STEROID_BL = 0).
    lcl <- log(0.184); label("Clearance CL for typical patient (L/day)")                     # Zheng 2016 Table 2, theta1
    lvc <- log(2.82);  label("Central volume of distribution V1 for typical patient (L)")    # Zheng 2016 Table 2, theta2
    lvp <- log(1.88);  label("Peripheral volume of distribution V2 (L)")                     # Zheng 2016 Table 2, theta3
    lq  <- log(0.34);  label("Inter-compartmental clearance Q (L/day)")                      # Zheng 2016 Table 2, theta4

    # Covariate effects on CL - power-form on WT and BGENE21, fractional on STEROID_BL.
    # CL = theta1 * (WT/64.3)^theta5 * (BGENE21/12.04)^theta6 * (1 + theta7 * STEROID_BL)
    e_wt_cl       <- 0.45; label("Power exponent of WT on CL (unitless, reference 64.3 kg)") # Zheng 2016 Table 2, theta5
    e_bgene21_cl  <- 0.09; label("Power exponent of BGENE21 on CL (unitless, reference 12.04)") # Zheng 2016 Table 2, theta6
    e_steroid_cl  <- 0.11; label("Fractional increase in CL for baseline steroid users")      # Zheng 2016 Table 2, theta7

    # Covariate effects on V1 - power-form on WT and DOSE, fractional on STEROID_BL.
    # V1 = theta2 * (WT/64.3)^theta8 * (DOSE/600)^theta9 * (1 - theta10 * STEROID_BL)
    e_wt_vc       <- 0.36; label("Power exponent of WT on V1 (unitless, reference 64.3 kg)") # Zheng 2016 Table 2, theta8
    e_dose_vc     <- 0.06; label("Power exponent of DOSE on V1 (unitless, reference 600 mg)") # Zheng 2016 Table 2, theta9
    e_steroid_vc  <- 0.09; label("Fractional decrease in V1 for baseline steroid users")      # Zheng 2016 Table 2, theta10

    # IIV - log-normal / exponential BSV assumed on CL and V1 (Zheng 2016 Table 2).
    # Converting CV% to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL BSV 24% -> omega^2 = log(0.24^2 + 1) = 0.056016
    #   V1 BSV 16% -> omega^2 = log(0.16^2 + 1) = 0.025278
    # No BSV estimated on Q or V2. The published final model also included inter-occasion
    # variability (IOV) on V1 of 26% CV for the first seven infusions; IOV is not carried in
    # this library implementation because occasion-specific randomness is outside the scope of
    # a typical single-eta subject-level model. See the vignette Assumptions section.
    etalcl ~ 0.056016  # Zheng 2016 Table 2, BSV CL = 24% CV
    etalvc ~ 0.025278  # Zheng 2016 Table 2, BSV V1 = 16% CV

    # Residual error - combined additive + proportional (Zheng 2016 Table 2).
    propSd <- 0.0701; label("Proportional residual error (fraction)")                        # Zheng 2016 Table 2, proportional error 7.01%
    addSd  <- 1.14;   label("Additive residual error (ug/mL)")                               # Zheng 2016 Table 2, additive error 1.14 ug/mL
  })

  model({
    # Covariate effects on CL and V1 per Zheng 2016 Equations (final model).
    cov_cl <- (WT / 64.3)^e_wt_cl * (BGENE21 / 12.04)^e_bgene21_cl * (1 + e_steroid_cl * STEROID_BL)
    cov_vc <- (WT / 64.3)^e_wt_vc * (DOSE / 600)^e_dose_vc         * (1 - e_steroid_vc * STEROID_BL)

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + etalvc) * cov_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    # Micro-constants for the two-compartment system with IV input into the central compartment.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system - IV dosing into the central compartment, first-order elimination from central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
