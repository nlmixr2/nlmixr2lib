Jaworowicz_2006_levalbuterol <- function() {
  description <- paste(
    "Two-compartment population PK model for (R)-albuterol following inhaled",
    "levalbuterol (90 ug) or racemic albuterol (180 ug) via a",
    "hydrofluoroalkane metered-dose inhaler in pediatric (4-11 years) and",
    "adult (12-81 years) asthma patients. First-order absorption, linear",
    "elimination, body-weight effects on apparent clearance (linear-additive)",
    "and central volume (power), and a pediatric-vs-adult split on absorption",
    "rate. The reference parameters are the Adult / Study 051-353 / single-dose",
    "levalbuterol-visit values (bioavailability anchor F1 = 1)."
  )
  reference <- paste(
    "Jaworowicz D, Maier G, Baumgartner RA, Hsu R, Grasela TH.",
    "Population pharmacokinetics of (R)-albuterol following inhaled",
    "levalbuterol or racemic albuterol via a hydrofluoroalkane metered",
    "dose inhaler in pediatric and adult asthma patients.",
    "Poster T3350, American Association of Pharmaceutical Scientists",
    "Annual Meeting and Exposition, San Antonio, TX, October 29 -",
    "November 2, 2006. Sepracor, Inc. / Cognigen Corporation."
  )
  vignette <- "Jaworowicz_2006_levalbuterol"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference body weight 74.8 kg (cohort median). Apparent clearance",
        "CL/F uses a linear-additive form CL = 59.1 + 0.477 * (WT - 74.8)",
        "L/hr (Equation 1 of poster). Apparent central volume Vc/F uses a",
        "power form Vc = 527 * (WT / 74.8)^0.361 L (Equation 2 of poster)."
      ),
      source_name        = "WTKG"
    ),
    CHILD = list(
      description        = "Pediatric age-cohort indicator (1 if subject < 12 years, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Used as the pediatric-vs-adult switch on absorption rate Ka:",
        "Adult subjects (>= 12 years) have Ka = 6.28 1/hr; pediatric",
        "subjects (< 12 years) have Ka = 3.08 1/hr. Encoded as a",
        "multiplicative log-ratio effect e_child_ka on lka."
      ),
      source_name        = "PED"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 632L,
    n_studies      = 3L,
    age_range      = "4-81 years (pediatric subset 4-11 years; adult subset 12-81 years)",
    age_median     = "Adults median ~36 years (433.7 +/- 192.4 months); pediatrics median ~9 years (106.9 +/- 28.1 months)",
    weight_range   = "14.5-167.5 kg",
    weight_median  = "75.2 kg overall (74.8 kg used as the reference for WT scaling); adults 80.8 kg, pediatrics 37.1 kg",
    sex_female_pct = 51,
    race_ethnicity = c(Caucasian = 69.6, Black = 19.3, Hispanic = 7.4, Asian = 2.4, Other = 1.3),
    disease_state  = "Adult and pediatric asthma",
    dose_range     = paste(
      "Inhaled HFA MDI: 90 ug levalbuterol QID or 180 ug racemic",
      "albuterol QID. Treatment duration 4 weeks (pediatrics) or 8 weeks",
      "(adults). Sampling at first dose (pre, 1-2, 4-6 hr) and after",
      "4 or 8 weeks (pre, 0.25, 0.5, 1, 2, 4, 8 hr)."
    ),
    regions        = "Three randomized multi-center placebo- and active-controlled Phase 3 trials (sponsor identifiers 051-353, 051-355, and a third pediatric trial); region not stated on the poster",
    notes          = paste(
      "PK measurand is (R)-albuterol plasma concentration (pg/mL in",
      "the source; this model emits ng/mL = 1e-3 pg/mL via Cc =",
      "central / vc with dose in ug and Vc in L). 429 subjects received",
      "levalbuterol HFA MDI and 203 received racemic albuterol HFA MDI;",
      "n = 3791 plasma samples. Baseline demographics in Table 1 of the",
      "poster. Bioavailability reference (F1 = 1) is the Study 051-353",
      "Visit 6 (single-dose levalbuterol) cohort; other strata reported",
      "as relative F1 values are documented in the vignette."
    )
  )

  ini({
    # Reference: Jaworowicz 2006 AAPS Poster T3350, Table 2 (Parameter
    # Estimates and Standard Errors for the Final PPK Model) and
    # Equations 1-2 (page footer of the Results column). Each value
    # below carries an in-file source-trace comment.

    # --- Structural absorption ---------------------------------------
    # Adult Ka is the reference; pediatric Ka enters as a log-ratio
    # covariate effect on lka via the CHILD indicator. log(3.08 / 6.28)
    # = -0.7126 (rounded to 4 d.p.).
    lka         <- log(6.28)
    label("Adult (>= 12 years) absorption rate constant ka (1/hr)")        # Table 2: Ka adults = 6.28 1/hr (%SEM 7.8)
    e_child_ka  <- log(3.08 / 6.28)
    label("Log-ratio effect of CHILD on lka (pediatric vs adult)")          # Table 2: Ka pediatrics = 3.08 1/hr (%SEM 15.4); log(3.08/6.28) = -0.7126

    # --- Apparent clearance: linear-additive WT scaling --------------
    # TVCL = 59.1 + 0.477 * (WT - 74.8) per Equation 1; encoded so the
    # log-normal IIV multiplies the entire TVCL (paper used NONMEM
    # EXP(ETA) on CL with TVCL constructed inside the structural model).
    lcl         <- log(59.1)
    label("Apparent clearance intercept at WT = 74.8 kg (L/hr)")             # Table 2: CL/F = 59.1 L/hr (%SEM 14.6); Equation 1 intercept
    e_wt_cl     <- 0.477
    label("Linear-additive WT effect on CL/F (L/hr per kg above 74.8 kg)")   # Table 2: slope = 0.477 L/hr/kg (%SEM 26.4); Equation 1 slope

    # --- Apparent central volume: power WT scaling -------------------
    lvc         <- log(527)
    label("Apparent central volume reference at WT = 74.8 kg (L)")           # Table 2: Vc/F = 527 L (%SEM 6.5); Equation 2 reference
    e_wt_vc     <- 0.361
    label("Power exponent of WT on Vc/F (unitless)")                         # Table 2: WT exponent = 0.361 (%SEM 27.2); Equation 2 exponent

    # --- Peripheral disposition --------------------------------------
    lvp         <- log(506)
    label("Apparent peripheral volume Vp (L)")                               # Table 2: Vp = 506 L (%SEM 25.1)
    lq          <- log(100)
    label("Inter-compartmental clearance Q (L/hr)")                          # Table 2: Q = 100 L/hr (%SEM 9.7); no IIV reported

    # --- Bioavailability reference -----------------------------------
    # F1 = 1.0 for the Adult / Study 051-353 / single-dose levalbuterol
    # visit (Table 2 footnote b). Other strata are documented in the
    # vignette: Adult LEV (Study 051-355) 0.707, Adult LEV (Study
    # 051-353 non-SD) 0.725, Pediatric LEV 0.550, Adult RAC (Study
    # 051-355) 1.01, Adult RAC (Study 051-353 SD) 0.880, Pediatric
    # RAC 0.830.
    lfdepot     <- fixed(log(1))
    label("Bioavailability anchor for depot (F1, reference cohort)")          # Table 2 footnote b: F1 = 1 reference (Adult LEV Study 051-353 SD visit)

    # --- IIV ---------------------------------------------------------
    # Reported as %CV in Table 2 ("Magnitude of IIV"). Convert to the
    # log-normal variance scale via omega^2 = log(CV^2 + 1). The
    # adult-cohort Ka IIV (67.60 %CV; omega^2 = 0.3764) is used as the
    # primary etalka variance; the pediatric Ka IIV (74.63 %CV;
    # omega^2 = 0.4424) is similar and documented in the vignette.
    etalka      ~ 0.3764                                                     # Table 2: Ka adults IIV 67.60 %CV -> omega^2 = log(1 + 0.676^2)
    etalcl      ~ 0.1450                                                     # Table 2: CL/F IIV 39.50 %CV -> omega^2 = log(1 + 0.395^2)
    etalvc      ~ 0.2078                                                     # Table 2: Vc/F IIV 48.06 %CV -> omega^2 = log(1 + 0.4806^2)
    etalvp      ~ 0.3927                                                     # Table 2: Vp   IIV 69.35 %CV -> omega^2 = log(1 + 0.6935^2)
    # F1 IIV reported per-stratum (Adult Study 051-355 LEV 24.31 %CV is
    # the only point estimate paired with an %SEM in Table 2; used here
    # as the etalfdepot variance). omega^2 = log(1 + 0.2431^2) = 0.0574.
    etalfdepot  ~ 0.0574                                                     # Table 2: F1 IIV 24.31 %CV (Adult Study 051-355 LEV) -> omega^2 = log(1 + 0.2431^2)

    # --- Residual error ----------------------------------------------
    # Table 2 reports the residual error as a NONMEM proportional-plus-
    # additive model in pg/mL with two components, "Additive RV (sigma2)"
    # and "Ratio of Additive/Proportional RV components (sigma2/sigma1)".
    # The footnote d states that the realised %CV is 21.99 % at IPRED =
    # 25 pg/mL and 21.33 % at IPRED = 900 pg/mL. Solving the standard
    # combined-error variance equation %CV^2 = sigma_prop^2 +
    # (sigma_add / IPRED)^2 from these two anchor points gives
    # sigma_prop = 0.2133 (proportional SD ~21.3 %) and sigma_add =
    # 1.34 pg/mL = 0.00134 ng/mL. The values entered below are this
    # reconstruction. See the vignette Assumptions section for the full
    # derivation; the raw Table 2 numbers (sigma2 = 0.0455, ratio
    # sigma2/sigma1 = 35.2) are NONMEM internal-scale variances and do
    # not transfer directly to the rxode2 SD-parameterised propSd /
    # addSd surface.
    propSd      <- 0.2133
    label("Proportional residual SD on Cc (fraction)")                       # Table 2 footnote d: proportional component from %CV(IPRED=900 pg/mL) = 21.33 %
    addSd       <- 0.00134
    label("Additive residual SD on Cc (ng/mL; 1.34 pg/mL in source units)")  # Table 2 footnote d: additive component derived from the %CV(25 pg/mL) - %CV(900 pg/mL) spread
  })

  model({
    # Individual PK parameters.
    # Ka splits adult vs pediatric via the CHILD indicator; the IIV is
    # shared (etalka) and uses the adult-cohort variance.
    ka          <- exp(lka + e_child_ka * CHILD + etalka)

    # Apparent clearance: linear-additive WT scaling around 74.8 kg.
    # IIV multiplies the WT-adjusted typical value (exponential-error
    # model with the structural TVCL inside the multiplication).
    cl          <- (exp(lcl) + e_wt_cl * (WT - 74.8)) * exp(etalcl)

    # Apparent central volume: power WT scaling around 74.8 kg.
    vc          <- exp(lvc + etalvc) * (WT / 74.8)^e_wt_vc

    # Peripheral volume and inter-compartmental clearance: no WT effect
    # reported; Q has no IIV.
    vp          <- exp(lvp + etalvp)
    q           <- exp(lq)

    # Micro-constants.
    kel         <- cl / vc
    k12         <- q  / vc
    k21         <- q  / vp

    # ODE system: depot -> central <-> peripheral1.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1)  <-                   k12     * central - k21 * peripheral1

    # Bioavailability anchor on the inhaled depot.
    f(depot)    <- exp(lfdepot + etalfdepot)

    # Observation: (R)-albuterol plasma concentration. Dose in ug, Vc
    # in L gives Cc in ug/L = ng/mL. The source poster reports Cc in
    # pg/mL; multiply Cc by 1000 to obtain pg/mL.
    Cc          <- central / vc

    Cc          ~ add(addSd) + prop(propSd)
  })
}
