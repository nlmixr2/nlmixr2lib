Li_2024_norepinephrine <- function() {
  description <- "Two-compartment population PK model with endogenous production and an input lag for norepinephrine in healthy volunteers, awake and under propofol/remifentanil general anesthesia (Li 2024)"
  reference <- "Li Y, Koomen JV, Eleveld DJ, van den Berg JP, Vos JJ, de Keijzer IN, Struys MMRF, Colin PJ. Population Pharmacokinetic Modelling of Norepinephrine in Healthy Volunteers Prior to and During General Anesthesia. Clin Pharmacokinet. 2024;63(11):1597-1608. doi:10.1007/s40262-024-01430-y"
  vignette <- "Li_2024_norepinephrine"

  # The model is parameterised in the units the paper reports its parameters
  # in: Prod in ng/min, CL and Q in L/min, Vc and Vp in L. The natural
  # concentration unit is therefore ng/L. The paper reports OBSERVED
  # concentrations in nmol/L; divide Cc by the norepinephrine molar mass
  # (169.18 ng/nmol) to reproduce the published nmol/L values.
  units <- list(time = "min", dosing = "ng", concentration = "ng/L")

  compartmentData <- list(
    central     = list(analyte = "norepinephrine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "norepinephrine", units = "ng", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling included a priori with theory-based fixed exponents (Methods 2.2): 0.75 on CL, Q and the endogenous production rate, 1.00 on Vc and Vp, 0.25 on the input lag time. Reference weight 70 kg.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL through Eq. 2, exp((theta_AGE~CL / 100) * (AGE - 35)); the centering age of 35 years is printed in the equation itself.",
      source_name        = "AGE"
    ),
    CONMED_PROPOFOL_CC = list(
      description        = "Measured plasma concentration of co-administered propofol",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = "0 (awake phase; no propofol administered)",
      notes              = "Time-varying covariate on CL (Results 3.2, Eq. 4). Set to 0 for the awake phase. In the general-anesthesia phase of the source study propofol was delivered by Eleveld-model target-controlled infusion to a 50% age-adjusted drug-effect concentration; the observed median measured concentration was 3.53 ug/mL. This covariate replaced the binary session factor F_SESS in the final model.",
      source_name        = "CPROP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36,
    n_studies      = 1,
    age_range      = "18-70 years",
    age_median     = "not reported overall; per-stratum medians 22, 39.5 and 57.5 years (mean 42, SD 16 years)",
    weight_range   = "51.5-94.8 kg",
    weight_median  = "not reported overall; per-stratum medians 66.9, 72.3 and 69.6 kg",
    sex_female_pct = 50,
    disease_state  = "healthy volunteers",
    dose_range     = "5 ug intravenous bolus followed by a five-step step-up infusion of 0.04, 0.08, 0.12, 0.16 and 0.20 ug/kg/min, each step 15 minutes; repeated in an awake phase and a general-anesthesia phase",
    regions        = "Netherlands (single centre, University Medical Center Groningen)",
    notes          = "Cross-over healthy-volunteer study (Netherlands Trial Register NL9312); 12 volunteers in each of three age strata (18-34, 35-50, 51-70 years) with a 1:1 male-to-female ratio. Baseline demographics are in Table 1. 1219 norepinephrine samples were analysed after excluding 4 outliers. Observed norepinephrine plasma concentrations ranged from 0.25 to 67.09 nmol/L. General anesthesia used propofol plus remifentanil by Eleveld-model target-controlled infusion; a 30-second 50 mA / 100 Hz tetanic electrical stimulus was applied at each dosing step as a surgical-incision surrogate."
  )

  ini({
    # Structural parameters (Table 3) -- typical values for a 70 kg, 35-year-old
    # volunteer with no propofol on board (CONMED_PROPOFOL_CC = 0).
    lcl   <- log(2.1)      ; label("Clearance (L/min)")                                    # Table 3, theta 1 (95% CI 2.073-2.201)
    lvc   <- log(2.4)      ; label("Central volume of distribution (L)")                   # Table 3, theta 2 (95% CI 1.988-2.850)
    lq    <- log(0.6)      ; label("Intercompartmental clearance (L/min)")                 # Table 3, theta 3 (95% CI 0.534-0.729); the table's "Q (min-1)" unit label is a typo -- Q is a clearance
    lvp   <- log(3.6)      ; label("Peripheral volume of distribution (L)")                # Table 3, theta 4 (95% CI 3.058-4.196)
    lkin  <- log(497.7)    ; label("Endogenous norepinephrine production rate (ng/min)")    # Table 3, theta 5 "Prod" (95% CI 487.905-545.814)
    ltlag <- log(13.7 / 60); label("Input lag time (min)")                                 # Table 3, theta 6 "Lag time" = 13.7 s; converted to minutes to match the model time unit

    # Allometric exponents -- theory-based and fixed a priori, not estimated
    # (Methods 2.2: "In accordance with theory based allometric scaling, the
    # exponents were 0.75 for clearance parameters and production rate, 1.00
    # for volumes of distribution parameters, and 0.25 for the lag time
    # parameter.").
    e_wt_cl   <- fixed(0.75); label("Allometric exponent on CL (unitless)")               # Methods 2.2
    e_wt_q    <- fixed(0.75); label("Allometric exponent on Q (unitless)")                # Methods 2.2
    e_wt_kin  <- fixed(0.75); label("Allometric exponent on the production rate (unitless)") # Methods 2.2
    e_wt_vc   <- fixed(1)   ; label("Allometric exponent on Vc (unitless)")               # Methods 2.2
    e_wt_vp   <- fixed(1)   ; label("Allometric exponent on Vp (unitless)")               # Methods 2.2
    e_wt_tlag <- fixed(0.25); label("Allometric exponent on the lag time (unitless)")     # Methods 2.2

    # Covariate effects on CL. Both coefficients are reported on a per-100
    # scale: Eq. 2 divides theta_AGE~CL by 100 explicitly, and applying the
    # same /100 to Eq. 4 reproduces the paper's own arithmetic exactly (see
    # the model() block and the vignette Errata).
    e_conmed_propofol_cc_cl <- -3.57 ; label("Effect of measured propofol concentration on CL (per 100 ug/mL, exponential)") # Table 3, theta 7 "CPROP~CL" (95% CI -4.312 to -2.780)
    e_age_cl                <- -0.344; label("Effect of age on CL (per 100 years, exponential)")                             # Table 3, theta 8 "AGE~CL" (95% CI -0.598 to -0.091)

    # IIV. Table 3 reports %CV; the table footnote gives the back-transform
    # omega^2 = log((CV/100)^2 + 1). IIV on Q was tested and removed
    # (Results 3.1, dOFV = 1.75).
    etalcl  ~ 0.0111733 ; label("IIV on clearance")                       # Table 3, IIV_CL   = 10.6% CV [shrinkage 9%];  log(0.106^2 + 1)
    etalvc  ~ 0.1806744 ; label("IIV on central volume")                  # Table 3, IIV_Vc   = 44.5% CV [shrinkage 13%]; log(0.445^2 + 1)
    etalvp  ~ 0.1511886 ; label("IIV on peripheral volume")               # Table 3, IIV_Vp   = 40.4% CV [shrinkage 12%]; log(0.404^2 + 1)
    etalkin ~ 0.0883918 ; label("IIV on the endogenous production rate")  # Table 3, IIV_prod = 30.4% CV [shrinkage 8%];  log(0.304^2 + 1)

    # Eq. 4 makes this eta ADDITIVE on the propofol coefficient,
    # (theta_CPROP~CL + eta_CPROP~CL), so it is a random effect on the
    # coefficient itself rather than a log-normal multiplier of it. Table 3
    # nonetheless reports it through the same %CV footnote formula, which is
    # what is inverted here.
    etae_conmed_propofol_cc_cl ~ 2.0901643; label("IIV on the propofol effect on clearance") # Table 3, IIV_CPROP~CL = 266.2% CV [shrinkage 32%]; log(2.662^2 + 1)

    # Combined residual error (Results 3.1). The footnote states that
    # within-subject variability is reported as sqrt(estimate), i.e. these
    # are already standard deviations.
    addSd  <- 32.4 ; label("Additive residual error (ng/L)")           # Table 3, dRUV additive = 32.4 (95% CI 27.2-38.1); ng/L, i.e. 0.19 nmol/L
    propSd <- 0.167; label("Proportional residual error (fraction)")   # Table 3, dRUV proportional = 16.7% (95% CI 15.7-17.6)
  })

  model({
    # Reference values used by the covariate models.
    wtRef  <- 70   # kg; standard allometric reference weight (Methods 2.2)
    ageRef <- 35   # years; centering age printed inside Eq. 2

    # Covariate effect on CL. Eq. 2: F_AGE = exp((theta_AGE~CL / 100) *
    # (AGE - 35)). Eq. 4: F_COV = exp((theta_COV~CL + eta_COV~CL) * (COV -
    # median baseline COV)); the baseline median of measured propofol is 0
    # because no propofol is present in the awake phase. Eq. 4 as printed
    # omits the /100 that Eq. 2 carries; the /100 is retained here because it
    # is the only reading that reproduces the paper's own stated effect sizes
    # (11.8% lower CL at the median propofol concentration of 3.53 ug/mL, and
    # ~3% lower CL per additional 10 years of age).
    fage  <- exp(e_age_cl / 100 * (AGE - ageRef))
    fprop <- exp((e_conmed_propofol_cc_cl + etae_conmed_propofol_cc_cl) / 100 * CONMED_PROPOFOL_CC)

    # Individual parameters
    cl   <- exp(lcl + etalcl)   * (WT / wtRef)^e_wt_cl * fage * fprop
    vc   <- exp(lvc + etalvc)   * (WT / wtRef)^e_wt_vc
    q    <- exp(lq)             * (WT / wtRef)^e_wt_q
    vp   <- exp(lvp + etalvp)   * (WT / wtRef)^e_wt_vp
    kin  <- exp(lkin + etalkin) * (WT / wtRef)^e_wt_kin
    tlag <- exp(ltlag)          * (WT / wtRef)^e_wt_tlag

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Endogenous baseline. Eq. 1 gives the central-compartment baseline
    # amount as Prod / (CL / Vc). The peripheral compartment is initialised
    # at the amount that makes the whole endogenous system stationary,
    # kin * vp / cl, which follows from setting both derivatives to zero;
    # without it the model produces a spurious redistribution dip in the
    # first ~20 minutes even though the paper's own baseline expression is
    # satisfied at t = 0.
    central(0)     <- kin * vc / cl
    peripheral1(0) <- kin * vp / cl

    d/dt(central)     <- kin - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Short delay between the infusion site and the contralateral-arm
    # sampling site (Methods 2.2; Results 3.1 restricted the estimate to the
    # 10-16 s circulation-time window).
    alag(central) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
