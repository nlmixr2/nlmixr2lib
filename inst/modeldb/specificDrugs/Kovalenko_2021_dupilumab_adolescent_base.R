Kovalenko_2021_dupilumab_adolescent_base <- function() {
  description <- "Dupilumab primary base population PK model for adolescents 12 to <18 years with moderate-to-severe atopic dermatitis (Kovalenko 2021): 2-compartment model parameterised in rates with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; body weight (reference 75 kg) is the only covariate, on central volume."
  reference <- "Kovalenko P, Kamal MA, Davis JD, Huniti N, Xu C, Bansal A, Shumel B, DiCioccio AT. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab in Adolescents and Children >=6 to <12 Years of Age Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2021;10(11):1345-1357. doi:10.1002/cpdd.986"
  vignette <- "Kovalenko_2021_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Primary BASE model for adolescents >=12 to <18 years (study R668-AD-1526,
  # NCT03054428), Kovalenko 2021 Table 1 and Supplementary Table 2
  # ('Adolescents >=12 to <18 Years of Age' column).  The structure is the
  # adult model of Kovalenko 2020 (doi:10.1002/cpdd.780) applied without
  # structural change (Methods; Figure 1).
  #
  # Estimation status (Methods): Vc, ke and Vm estimated; kcp, kpc, ka,
  # MTT and F fixed to the adult values; Km fixed at 0.01 mg/L (assessed by
  # log-likelihood profiling as in the adults).
  #
  # Transit chain: Figure 1 shows the SC dose entering a chain of three
  # transit compartments spanned by the mean transit time (MTT), which feed
  # the absorption depot drained at ka.  The three MTT boxes are depot,
  # transit1 and transit2 here, and the ka depot is transit3; ktr = 3 / MTT
  # so that MTT is the mean time spent in the chain (the NN = 3,
  # KTR = NN / MTT convention of the same group's later deposited dupilumab
  # control stream, Nguyen 2026 doi:10.1002/cpt.70233).
  compartmentData <- list(
    depot = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "dupilumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dupilumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "dupilumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on central volume Vc, (WT/75)^e_wt_vc.  Methods: 'When weight was explored as a covariate, the central value was set to 75 kg in the primary analyses' (chosen for comparability across the age populations, not the cohort median).",
      source_name = "weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 162L,
    n_studies = 1L,
    age_range = "12 to <18 years",
    age_mean = "14.4 (SD 1.59) years",
    weight_mean = "65.3 (SD 22.0) kg",
    sex_female_pct = 43.3,
    disease_state = "Moderate-to-severe atopic dermatitis, dupilumab monotherapy (LIBERTY AD ADOL, R668-AD-1526, NCT03054428).",
    dose_range = "SC dupilumab 200 mg q2w (<60 kg, n = 43) or 300 mg q2w (>=60 kg, n = 39), or 300 mg q4w irrespective of weight (n = 82).",
    regions = "United States and Canada (Supplementary Table 1).",
    notes = "162 of 165 adolescents on active treatment and 827 of 1006 samples in the primary analysis (Methods); 56.7% male.  Sparse sampling on days 1, 29, 57, 113 and end of treatment.  BLQ data handled by Beal M3.  Estimation in Monolix 2019R2 (SAEM + importance sampling)."
  )

  ini({
    # Estimated structural parameters (Table 1 / Supp. Table 2, adolescent column)
    lvc   <- log(2.54);   label("Central volume at 75 kg (L)")                                   # Table 1 adolescents: Vc = 2.54 L (SE 0.0473)
    lkel  <- log(0.0508); label("Linear elimination rate constant ke (1/d)")                     # Table 1 adolescents: ke = 0.0508 1/d (SE 0.00172)
    lvmax <- log(1.46);   label("Maximum target-mediated elimination rate Vm (mg/L/d)")          # Table 1 adolescents: Vm = 1.46 mg/L/d (SE 0.0314)

    # Fixed structural parameters
    km      <- fixed(0.01);       label("Michaelis-Menten constant Km (mg/L)")                  # Table 1 adolescents: Km = 0.01 (fixed)
    lkcp    <- fixed(log(0.211)); label("Central-to-peripheral rate constant kcp (1/d)")        # Table 1 adolescents: kcp = 0.211 (fixed)
    lkpc    <- fixed(log(0.310)); label("Peripheral-to-central rate constant kpc (1/d)")        # Table 1 adolescents: kpc = 0.310 (fixed)
    lka     <- fixed(log(0.306)); label("Absorption rate constant ka (1/d)")                    # Table 1 adolescents: ka = 0.306 (fixed)
    lmtt    <- fixed(log(0.105)); label("Mean transit time MTT (d)")                            # Table 1 adolescents: MTT = 0.105 (fixed)
    lfdepot <- fixed(log(0.642)); label("Subcutaneous bioavailability F (fraction)")            # Table 1 adolescents: F = 0.642 (fixed)

    # Covariate effect
    e_wt_vc <- 0.853; label("Power exponent of WT/75 on Vc (unitless)")                          # Table 1 adolescents: Vc ~ weight = 0.853 (SE 0.0438)

    # IIV (Supp. Table 2 adolescents): SD(ln Vc) = 0.141, SD(ln ke) = 0.335,
    # Corr(ln ke, ln Vc) = -0.407.  Variances 0.141^2, 0.335^2;
    # covariance -0.407 * 0.141 * 0.335.
    etalvc + etalkel ~ c(0.019881,
                         -0.0192246, 0.112225)

    # Residual error (Supp. Table 2 adolescents)
    propSd <- 0.0990; label("Proportional residual error (fraction)")                            # Supp. Table 2 adolescents: sigma_prop = 9.90 CV% (SE 0.593)
    addSd  <- 2.41;   label("Additive residual error (mg/L)")                                    # Supp. Table 2 adolescents: sigma_add = 2.41 mg/L (SE 0.248)
  })
  model({
    vc   <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    kel  <- exp(lkel + etalkel)
    vmax <- exp(lvmax)
    kcp  <- exp(lkcp)
    kpc  <- exp(lkpc)
    ka   <- exp(lka)
    mtt  <- exp(lmtt)

    # Three transit transfers at ktr (depot -> transit1 -> transit2 -> transit3);
    # transit3 is the Figure 1 absorption depot drained at ka.
    ktr <- 3 / mtt

    Cc <- central / vc

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * (depot - transit1)
    d/dt(transit2)    <-  ktr * (transit1 - transit2)
    d/dt(transit3)    <-  ktr * transit2 - ka * transit3
    d/dt(central)     <-  ka * transit3 - kel * central - kcp * central + kpc * peripheral1 -
                          vc * vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1

    # F applies to the SC dose; an IV dose into central bypasses it
    f(depot) <- exp(lfdepot)

    Cc ~ add(addSd) + prop(propSd)
  })
}
