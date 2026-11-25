Singh_2016_stamulumab <- function() {
  description <- "Two-compartment population PK model with linear elimination for stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, in healthy adult volunteers after single ascending IV doses of 0.3-11.25 mg/kg (Singh 2016). Parameters are reported per kg body weight and scale linearly with WT; log-normal IIV on CL and V1 and proportional residual error. The paper uses it to project steady-state exposure at the 1, 3 and 10 mg/kg every-2-week doses of the phase I/II muscular dystrophy trial."
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
    species = "human",
    n_subjects = 61,
    n_studies = 1,
    age_range = "18-45 years",
    weight_range = NA,
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state = "Healthy adult volunteers (US phase I single ascending dose study)",
    dose_range = "Single IV doses of 0.3, 1, 3, 7.5 and 11.25 mg/kg, or placebo",
    regions = "United States",
    notes = "Singh 2016 Methods 'MYO-029 clinical trials': randomized, double-blind, placebo-controlled, sequential-group SAD study (n = 61 including placebo). PK sampled predose and at 0.5, 1, 4, 12, 24, 48, 72, 168, 312, 480, 648, 816, 984, 1152 and 1320 h after the start of the infusion. The infusion duration of the SAD study is not reported; the paper's multiple-dose projections (Supplementary Table S1) use a 1-h infusion. Serum MYO-029 by GDF-8-capture ELISA (63.2-720 ng/mL in 100% human serum)."
  )

  ini({
    # Singh 2016 Table 1, 'Human (%RSE)' column; RSE in parentheses.
    lcl <- log(0.38); label("Clearance per kg body weight (mL/h/kg)") # Table 1, CL 'Clearance from central compartment' human = 0.38 mL/h/kg (4% RSE)
    lvc <- log(65); label("Central volume per kg body weight (mL/kg)") # Table 1, V1 'Volume of central compartment' human = 65 mL/kg (3% RSE)
    lq <- log(0.31); label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table 1, Q 'Distributive clearance' human = 0.31 mL/h/kg (12% RSE)
    lvp <- log(33); label("Peripheral volume per kg body weight (mL/kg)") # Table 1, V2 'Volume of peripheral compartment' human = 33 mL/kg (6% RSE)

    # Table 1 prints the IIV as omega^2 (variances of log-normal etas).
    etalcl ~ 0.04 # Table 1, 'omega CL(IIV)^2' human = 0.04 (27% RSE)
    etalvc ~ 0.03 # Table 1, 'omega V1(IIV)^2' human = 0.03 (25% RSE)

    # Table 1 prints sigma1 = 0.16 (not squared, unlike the omegas). Read as
    # the proportional SD: with it the model reproduces the 5th-95th
    # percentile of the peak in Supplementary Table S1, whereas a variance of
    # 0.16 (SD 0.4) would roughly double that interval (see the vignette).
    propSd <- 0.16; label("Proportional residual error (fraction)") # Table 1, sigma1 'Proportional error' human = 0.16 (8% RSE)
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
