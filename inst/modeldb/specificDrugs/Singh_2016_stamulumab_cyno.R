Singh_2016_stamulumab_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model with linear elimination for stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, after weekly IV doses of 10, 30 and 100 mg/kg in the 5-week and 39-week toxicology studies (Singh 2016). Parameters are reported per kg body weight and scale linearly with WT; log-normal IIV on CL and V1. Residual error is not reported."
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
    species = "cynomolgus monkey",
    n_subjects = NA,
    n_studies = 2,
    age_range = NA,
    weight_range = NA,
    sex_female_pct = 50,
    race_ethnicity = NA,
    disease_state = "Healthy male and female cynomolgus monkeys (toxicology studies)",
    dose_range = "10, 30 and 100 mg/kg IV once weekly for 5 weeks (5-week study) or 39 weeks (39-week study)",
    regions = NA,
    notes = "Singh 2016 Methods: 5-week study sampled after the day 1 and day 22 doses (10 min, 1, 8, 24, 48, 72 and 168 h) plus a day 15 trough; 39-week study (n = 4/sex/group) sampled over the day 169 dosing interval (0.083, 6, 24, 48, 72, 120 and 168 h). The 5-week study group size is not reported, so the total number of animals is unknown. Serum MYO-029 by GDF-8-capture ELISA (63.2-720 ng/mL in 100% monkey serum)."
  )

  ini({
    # Singh 2016 Table 1, 'Monkey (%RSE)' column; RSE in parentheses.
    lcl <- log(0.228); label("Clearance per kg body weight (mL/h/kg)") # Table 1, CL 'Clearance from central compartment' monkey = 0.228 mL/h/kg (4% RSE)
    lvc <- log(41.8); label("Central volume per kg body weight (mL/kg)") # Table 1, V1 'Volume of central compartment' monkey = 41.8 mL/kg (6% RSE)
    lq <- log(0.815); label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table 1, Q 'Distributive clearance' monkey = 0.815 mL/h/kg (6% RSE)
    lvp <- log(27.8); label("Peripheral volume per kg body weight (mL/kg)") # Table 1, V2 'Volume of peripheral compartment' monkey = 27.8 mL/kg (11% RSE)

    # Table 1 prints the IIV as omega^2 (variances of log-normal etas).
    etalcl ~ 0.27 # Table 1, 'omega CL(IIV)^2' monkey = 0.27 (2% RSE)
    etalvc ~ 0.35 # Table 1, 'omega V1(IIV)^2' monkey = 0.35 (3% RSE)

    # Residual error is not reported for the monkey fit (Table 1 '-'); fixed to zero.
    propSd <- fixed(0); label("Proportional residual error (fraction)") # not reported in Table 1 (monkey column '-')
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
