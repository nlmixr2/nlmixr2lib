Singh_2016_stamulumab_rat <- function() {
  description <- "Preclinical (rat, Sprague-Dawley). Two-compartment population PK model with linear elimination for stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, after single 15-minute IV infusions of 2, 10 and 50 mg/kg in male Sprague-Dawley rats (Singh 2016). Parameters are reported per kg body weight and scale linearly with WT; log-normal IIV on CL and V1. Residual error is not reported."
  reference <- "Singh P, Rong H, Gordi T, Bosley J, Bhattacharya I. Translational Pharmacokinetic/Pharmacodynamic Analysis of MYO-029 Antibody for Muscular Dystrophy. Clin Transl Sci. 2016;9(6):302-310. doi:10.1111/cts.12420. PMID 27700008. PMCID PMC5351001. MYO-029 is the development code of the antibody assigned the INN stamulumab."
  vignette <- "Singh_2016_stamulumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Singh 2016 Table 1 reports every parameter per kg body weight",
        "(mL/h/kg, mL/kg), so CL, V1, Q and V2 scale linearly with WT and",
        "doses are given as mg/kg * WT. No reference weight is needed."
      ),
      source_name = "per-kg normalisation (no WT column printed)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "stamulumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "stamulumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 18,
    n_studies = 1,
    age_range = NA,
    weight_range = NA,
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state = "Healthy male Sprague-Dawley rats",
    dose_range = "Single 15-min IV infusion of 2, 10 or 50 mg/kg (n = 6 per dose)",
    regions = NA,
    notes = "Singh 2016 Methods 'Rat pharmacokinetics studies': serum sampled predose and at 0.25, 1, 4, 24, 48, 72, 168, 336, 504, 672, 840, 1008 and 1056 h. Serum MYO-029 by GDF-8-capture ELISA (21.1-240 ng/mL in 100% rat serum)."
  )

  ini({
    # Singh 2016 Table 1, 'Rat (%RSE)' column; RSE in parentheses.
    lcl <- log(0.542); label("Clearance per kg body weight (mL/h/kg)") # Table 1, CL 'Clearance from central compartment' rat = 0.542 mL/h/kg (4% RSE)
    lvc <- log(58.9); label("Central volume per kg body weight (mL/kg)") # Table 1, V1 'Volume of central compartment' rat = 58.9 mL/kg (4% RSE)
    lq <- log(1.79); label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table 1, Q 'Distributive clearance' rat = 1.79 mL/h/kg (9% RSE)
    lvp <- log(95.3); label("Peripheral volume per kg body weight (mL/kg)") # Table 1, V2 'Volume of peripheral compartment' rat = 95.3 mL/kg (5% RSE)

    # Table 1 prints the IIV as omega^2 (variances of log-normal etas).
    etalcl ~ 0.15 # Table 1, 'omega CL(IIV)^2' rat = 0.15 (1% RSE)
    etalvc ~ 0.12 # Table 1, 'omega V1(IIV)^2' rat = 0.12 (1% RSE)

    # Residual error is not reported for the rat fit (Table 1 '-'); fixed to zero.
    propSd <- fixed(0); label("Proportional residual error (fraction)") # not reported in Table 1 (rat column '-')
  })

  model({
    # Per-kg parameters scaled to the individual: mL -> L.
    cl <- exp(lcl + etalcl) * WT / 1000
    vc <- exp(lvc + etalvc) * WT / 1000
    q <- exp(lq) * WT / 1000
    vp <- exp(lvp) * WT / 1000

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
