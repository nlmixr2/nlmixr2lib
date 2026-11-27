Obert_2026_fentanyl_rat <- function() {
  description <- "Preclinical (rat). Population PK model for intravenous fentanyl and its metabolite norfentanyl in Sprague-Dawley rats given repeated 100 ug/kg IV boluses twice weekly for 5 weeks: two-compartment fentanyl with linear clearance, Michaelis-Menten (non-depleting) norfentanyl formation from the fentanyl central concentration into a one-compartment norfentanyl disposition model; female sex lowers fentanyl clearance, norfentanyl formation Vmax and norfentanyl clearance, and fentanyl clearance is lower from the 10th injection onward."
  reference <- paste(
    "Obert DP, McKnite A, Park GH, Watt KM, Reilly CA, Vincent KF, Rower JE, Solt K.",
    "Pharmacokinetics of repeated, high-dose fentanyl administration in rats.",
    "Drug Metab Dispos. 2026;54(9):100390. doi:10.1016/j.dmd.2026.100390"
  )
  vignette <- "Obert_2026_fentanyl_rat"
  units <- list(time = "h", dosing = "ng", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "The control stream (Supplemental Material 1) codes SEX = 1 as the",
        "reference ('Most common') and applies (1 + THETA) when SEX = 0; the",
        "paper text states fentanyl clearance, Vmax and norfentanyl clearance",
        "are lower in females, so SEX = 0 is female and SEXF = 1 - SEX.",
        "The negative coefficients therefore carry over unchanged."
      ),
      source_name = "SEX"
    ),
    CYCLE = list(
      description = "Injection (dose) number of the twice-weekly fentanyl course, 1 = first injection",
      units = "(count)",
      type = "count",
      reference_category = "injections 1-9 (CYCLE < 10)",
      notes = paste(
        "Source column INJ, which took the values 1, 5 and 10 (the sampled",
        "injections; 4 animals were sampled at the 11th injection instead of",
        "the 10th, and the control stream handles only INJ = 1, 5 and 10, so",
        "those samples were coded 10). Increment CYCLE at each injection.",
        "The model dichotomizes it as in the final control stream:",
        "cycle10 = (CYCLE >= 10) multiplies fentanyl clearance by",
        "(1 + e_cycle10_cl). Only injections 1, 5 and 10 informed the estimate,",
        "so the factor applied to injections 2-4 and 6-9 (none) is an",
        "extrapolation."
      ),
      source_name = "INJ"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (as grams, WTG) on the base model but not retained. Body weight still sets the 100 ug/kg dose amount."
    ),
    AGE = list(
      description = "Age",
      units = "weeks",
      type = "continuous",
      notes = "Tested in the stepwise covariate search but not retained."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "fentanyl",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "fentanyl",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    ),
    central_norfentanyl = list(
      analyte = "norfentanyl",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "adult; mean baseline age 13.4 weeks (females) and 8.9 weeks (males)",
    weight_range = "280-380 g at purchase",
    weight_median = "mean baseline weight 286 g (females) and 343 g (males)",
    sex_female_pct = 47.4,
    disease_state = "healthy, opioid-naive at the first injection; repeated overdose-level fentanyl exposure",
    dose_range = "100 ug/kg IV bolus (tail vein) twice weekly, at least 72 h apart, for 5 weeks (10 injections)",
    regions = "United States (Massachusetts General Hospital)",
    notes = paste(
      "9 female and 10 male rats, giving 7 animals per sex with data at each",
      "sampled injection (Results, Study population). Arterial plasma sampled",
      "at 3, 15, 30, 60, 180 and 360 min after the 1st, 5th and 10th",
      "injections (11th for 4 animals). Fentanyl and norfentanyl assayed by",
      "UHPLC-ESI-MS/MS, calibrated range 0.01-10 ng/mL with dilution above."
    )
  )

  ini({
    # Fentanyl disposition (Table 3; male, injection 1 or 5 reference)
    lcl <- log(774); label("Fentanyl clearance, CLFent (mL/h)") # Table 3 'CLFent (mL/h)' = 774 (5% RSE)
    lvc <- log(275); label("Fentanyl central volume, VFent,Central (mL)") # Table 3 'VFent,Central (mL)' = 275 (9% RSE)
    lq <- log(232); label("Fentanyl intercompartmental clearance, QFent (mL/h)") # Table 3 'QFent (mL/h)' = 232 (18% RSE)
    lvp <- log(326); label("Fentanyl peripheral volume, VFent,Peripheral (mL)") # Table 3 'VFent,Peripheral (mL)' = 326 (9% RSE)

    # Norfentanyl formation (Michaelis-Menten) and disposition
    lkm_norfentanyl <- log(2.89); label("Michaelis constant of norfentanyl formation, on the fentanyl central concentration (ng/mL)") # Table 3 'Km (ng/mL)' = 2.89 (17% RSE)
    lvmax_norfentanyl <- log(78.8); label("Maximum rate of norfentanyl formation, Vmax (ng fentanyl-equivalents/h)") # Table 3 'Vmax' = 78.8 (12% RSE); printed as ng/mL per h but the control stream adds it to dA3/dt, an amount rate
    lcl_norfentanyl <- log(23.6); label("Norfentanyl clearance, CLNor (mL/h)") # Table 3 'CLNor (mL/h)' = 23.6 (11% RSE)
    lvc_norfentanyl <- log(10.4); label("Norfentanyl volume, VNor (mL)") # Table 3 'VNor (mL)' = 10.4 (18% RSE)

    # Covariate effects, all of the form (1 + theta * indicator)
    e_sexf_cl <- -0.158; label("Fractional change in fentanyl clearance in females (unitless)") # Table 3 'CLFent:Sex' = -0.158 (26% RSE)
    e_cycle10_cl <- -0.0796; label("Fractional change in fentanyl clearance from the 10th injection (unitless)") # Table 3 'CLFent:Injection' = -0.0796 (54% RSE)
    e_sexf_vmax_norfentanyl <- -0.764; label("Fractional change in norfentanyl formation Vmax in females (unitless)") # Table 3 'Vmax:Sex' = -0.764 (7% RSE)
    e_sexf_cl_norfentanyl <- -0.445; label("Fractional change in norfentanyl clearance in females (unitless)") # Table 3 'CLNor:Sex' = -0.445 (27% RSE)

    # IIV (exponential, variances; Table 3)
    etalcl ~ 0.00843 # Table 3 'omega2 CL,Fent' = 0.00843 (20% RSE, 9% shrinkage)
    etalvmax_norfentanyl ~ 0.0716 # Table 3 'omega2 Vmax' = 0.0716 (37% RSE, 11% shrinkage)
    etalcl_norfentanyl ~ 0.0602 # Table 3 'omega2 CL,Nor' = 0.0602 (65% RSE, 24% shrinkage)

    # Residual error: NONMEM $SIGMA variances of Y = F*(1 + ERR), SD = sqrt(variance)
    propSd <- 0.33466; label("Proportional residual error, fentanyl (fraction)") # Table 3 'Prop. error, fentanyl' = 0.112 (variance); sqrt(0.112) = 0.33466
    propSd_norfentanyl <- 0.51088; label("Proportional residual error, norfentanyl (fraction)") # Table 3 'Prop. error, norfentanyl' = 0.261 (variance); sqrt(0.261) = 0.51088
  })

  model({
    # Injection-number indicator, dichotomized as in the final control stream
    # (INJ 10 vs INJ 1 or 5)
    cycle10 <- CYCLE >= 10

    cl <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF) * (1 + e_cycle10_cl * cycle10)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    vmax_norfentanyl <- exp(lvmax_norfentanyl + etalvmax_norfentanyl) * (1 + e_sexf_vmax_norfentanyl * SEXF)
    km_norfentanyl <- exp(lkm_norfentanyl)
    cl_norfentanyl <- exp(lcl_norfentanyl + etalcl_norfentanyl) * (1 + e_sexf_cl_norfentanyl * SEXF)
    vc_norfentanyl <- exp(lvc_norfentanyl)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_norfentanyl <- cl_norfentanyl / vc_norfentanyl

    Cfent <- central / vc
    # Formation is driven by, but not subtracted from, the fentanyl central
    # compartment: DADT(1) in the control stream carries only CL and Q terms.
    form_norfentanyl <- vmax_norfentanyl * Cfent / (km_norfentanyl + Cfent)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_norfentanyl) <- form_norfentanyl - kel_norfentanyl * central_norfentanyl

    # central_norfentanyl is in fentanyl-mass equivalents (the dose is fentanyl
    # ng); the control stream scales it by S3 = V3 * (336.5 / 232.3), the
    # fentanyl:norfentanyl molecular-weight ratio, to report norfentanyl ng/mL.
    Cc <- Cfent
    Cc_norfentanyl <- central_norfentanyl / vc_norfentanyl * (232.3 / 336.5)

    Cc ~ prop(propSd)
    Cc_norfentanyl ~ prop(propSd_norfentanyl)
  })
}
