Williams_2012_valproic_acid_pediatric <- function() {
  description <- "Two-compartment population PK model for valproic acid in pediatric patients with epilepsy (Williams 2012). Allometric weight scaling on CL/Q (fixed 0.75) and Vc/Vp (fixed 1.0); estimated age power (-0.267) on Vc; reference weight 70 kg, reference age 8.5 years. Default first-order oral absorption is for divalproex sodium enteric-coated sprinkle (Ka 1.2 1/h, ALAG 1 h, FIXED); other formulations (syrup K0=410 mg/h, capsule Ka=2 1/h, tablet Ka=4.1 1/h with ALAG=2 h) require overriding lka/ltlag at simulation time. Direct IV dosing into the central compartment is supported. Residual error defaults to the TDM-subset proportional SD (CV 34.8%); paper also reports a TRIAL-subset SD (CV 4.6%) for the IV-infusion subset."
  reference <- "Williams JH, Jayaraman B, Swoboda KJ, Barrett JS. Population pharmacokinetics of valproic acid in pediatric patients with epilepsy: considerations for dosing spinal muscular atrophy patients. J Clin Pharmacol. 2012;52(11):1676-1688. doi:10.1177/0091270011428138. PMID 22167565."
  vignette <- "Williams_2012_valproic_acid_pediatric"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference weight 70 kg. Fixed exponents 0.75 on CL/Q and 1.0 on Vc/Vp (Williams 2012 Table I).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law effect on Vc with reference 8.5 years (median age in epilepsy data set). Estimated exponent -0.267 (Williams 2012 Table I, Final Model).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52,
    n_studies      = 2,
    age_range      = "1-17 years",
    age_median     = "8.5 years",
    weight_range   = "not explicitly tabulated; median 27.5 kg",
    weight_median  = "27.5 kg",
    sex_female_pct = 30.8,
    race_ethnicity = "Not retained as covariate (>10% missing in the analysis dataset).",
    disease_state  = "Pediatric epilepsy (10 subjects from an IV-infusion clinical trial of Depacon; 42 subjects from therapeutic drug monitoring at The Children's Hospital of Philadelphia, 2004-2006).",
    dose_range     = "Clinical-trial IV infusion 14 mg/kg (range 12-15 mg/kg) single dose; TDM-subset oral or IV multiple doses 23 mg/kg/day (range 3-60 mg/kg/day) across syrup, capsule, divalproex sprinkle, and tablet formulations.",
    regions        = "United States (multicenter IV infusion trial; Children's Hospital of Philadelphia TDM cohort).",
    notes          = "Final analysis data set: 231 observations across 52 subjects with 1-15 observations per subject. 36 male / 16 female. 13 subjects classified as induced (concomitant antiepileptic drugs); 39 monotherapy. Extended-release tablet dosing and post-15-hour TDM observations excluded. See Williams 2012 Results section 'Epilepsy Patient Population and Data Characteristics' and Table I."
  )

  ini({
    # ----------------------------------------------------------------
    # Structural parameters (Williams 2012 Table I, Final Model)
    # Reference subject: 70 kg, 8.5 years
    # ----------------------------------------------------------------
    lcl  <- log(0.854); label("Clearance for a 70 kg subject (L/h)")                  # Williams 2012 Table I (CL = THETA7 = 0.854 L/h, 6.21% RSE)
    lvc  <- log(10.3);  label("Central volume for a 70 kg, 8.5 yr subject (L)")       # Williams 2012 Table I (Vc = THETA8 = 10.3 L, 6.77% RSE)
    lq   <- log(5.34);  label("Intercompartmental clearance for a 70 kg subject (L/h)") # Williams 2012 Table I (Q = THETA9 = 5.34 L/h, 20.0% RSE)
    lvp  <- log(4.08);  label("Peripheral volume for a 70 kg subject (L)")            # Williams 2012 Table I (Vp = THETA10 = 4.08 L, 17.2% RSE)

    # ----------------------------------------------------------------
    # Absorption parameters - FIXED to literature values per Williams
    # 2012 Methods/Results: "absorption parameters [were fixed] to
    # previously characterized rates ... shown in Table I." Defaults
    # below are for divalproex sodium enteric-coated sprinkle.
    # Alternative formulations (override lka/ltlag at simulation time):
    #   Syrup:   K0 = 410 mg/h (zero-order; set dose rate accordingly)
    #   Capsule: Ka = 2.0 1/h, no lag
    #   Tablet:  Ka = 4.1 1/h, ALAG = 2 h
    # All four oral formulations approximate 100% bioavailability.
    # ----------------------------------------------------------------
    lka     <- fixed(log(1.2));   label("First-order absorption rate constant for sprinkle (1/h, FIXED)") # Williams 2012 Table I (Sprinkle Ka = THETA3 = 1.2 1/h)
    ltlag   <- fixed(log(1));     label("Absorption lag time for sprinkle (h, FIXED)")                    # Williams 2012 Table I (Sprinkle ALAG = THETA4 = 1 h)
    lfdepot <- fixed(log(1));     label("Oral bioavailability (fraction, FIXED at 1)")                    # Williams 2012 Discussion: "approximately 100% bioavailability" for all 4 oral formulations

    # ----------------------------------------------------------------
    # Allometric weight exponents (Williams 2012 Table I, theoretical
    # allometric values; reported without RSE in the Final Model column,
    # consistent with FIXED estimation)
    # ----------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless, FIXED)")  # Williams 2012 Table I (WT power CL = 0.75, WT power Q = 0.75)
    e_wt_vc_vp <- fixed(1.0);  label("Shared allometric exponent on Vc and Vp (unitless, FIXED)") # Williams 2012 Table I (WT power Vc = 1.0, WT power Vp = 1.0)

    # ----------------------------------------------------------------
    # Age effect on Vc (estimated; bootstrap 95% CI -0.378 to 0.211
    # spans zero per Williams 2012, "interpret with caution")
    # ----------------------------------------------------------------
    e_age_vc <- -0.267; label("Power-law age exponent on Vc relative to AGE=8.5 (unitless)")      # Williams 2012 Table I (AGE power Vc = THETA11 = -0.267, 18.2% RSE)

    # ----------------------------------------------------------------
    # IIV - Full block OMEGA on CL, Vc, Vp (Williams 2012 Table I, Final
    # Model). Variances on the lower triangle, row-major:
    #   omega^2_11 = var(CL),    omega_21 = cov(CL,Vc), omega^2_22 = var(Vc),
    #   omega_31 = cov(CL,Vp),   omega_32 = cov(Vc,Vp), omega^2_33 = var(Vp).
    # Bootstrap CIs for the covariance terms span zero -- Williams
    # 2012 cautions the off-diagonals are imprecisely estimated.
    # ----------------------------------------------------------------
    etalcl + etalvc + etalvp ~ c(0.129,
                                  0.0397, 0.0384,
                                  0.0777, 0.144,  1.02)
    # Williams 2012 Table I (omega^2_11 = 0.129 CV%35.9; omega_21 = 0.0397;
    # omega^2_22 = 0.0384 CV%19.6; omega_31 = 0.0777; omega_32 = 0.144;
    # omega^2_33 = 1.02 CV%101.5).

    # ----------------------------------------------------------------
    # Residual error - proportional. Williams 2012 reported separate
    # values for the TDM subset (sigma^2_prop = 0.121, CV% 34.8) and the
    # IV-infusion clinical-trial subset (sigma^2_prop = 0.00214, CV% 4.6).
    # The model file defaults to the TDM value (the dominant data
    # source, 174 / 231 observations). The TRIAL value is documented
    # here for users replicating Figure 2A or simulating prospective
    # IV clinical-trial scenarios.
    # ----------------------------------------------------------------
    propSd <- 0.348; label("Proportional residual SD, TDM subset (fraction; CV 34.8%)")           # Williams 2012 Table I (sigma^2_prop[TDM] = 0.121, sqrt = 0.348)
    # TRIAL-subset propSd = sqrt(0.00214) = 0.0463 (CV 4.6%)
  })

  model({
    # 1. Individual PK parameters with allometric WT scaling and AGE
    #    power on Vc (Williams 2012 Final Model, Table I)
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * (AGE / 8.5)^e_age_vc
    q   <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp  <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    tlag <- exp(ltlag)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. Two-compartment ODE system with first-order oral absorption
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Bioavailability and lag-time for the oral depot
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # 5. Observation and residual error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
