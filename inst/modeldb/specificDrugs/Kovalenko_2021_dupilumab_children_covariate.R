Kovalenko_2021_dupilumab_children_covariate <- function() {
  description <- "Dupilumab primary covariate population PK model for children 6 to <12 years with severe atopic dermatitis (Kovalenko 2021): 2-compartment model parameterised in rates with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; body weight and serum albumin on central volume and EASI score on the linear elimination rate."
  reference <- "Kovalenko P, Kamal MA, Davis JD, Huniti N, Xu C, Bansal A, Shumel B, DiCioccio AT. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab in Adolescents and Children >=6 to <12 Years of Age Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2021;10(11):1345-1357. doi:10.1002/cpdd.986. EASI centring value from the companion trial's baseline table: Siegfried EC, et al. Am J Clin Dermatol. 2023;24(5):787. doi:10.1007/s40257-023-00791-7 (post hoc analysis of LIBERTY AD PEDS), Table 1."
  vignette <- "Kovalenko_2021_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Primary COVARIATE model for children >=6 to <12 years (study
  # R668-AD-1652, NCT03345914), Kovalenko 2021 Table 2 and Supplementary
  # Table 3 ('Children >=6 to <12 Years of Age' column).  Same structure as
  # Kovalenko_2021_dupilumab_children_base.R.
  #
  # Estimation status (Methods): Vc and ke estimated; Vm = 1.64 mg/L/d was
  # estimated in the base model and FIXED in the covariate model (Table 2
  # footnote: a sensitivity re-estimation gave 1.57 mg/L/d); ka = 0.641 1/d
  # fixed from the phase 2a R668-AD-1412 fit; kcp, kpc, MTT, F, Km fixed to
  # the adult values.
  #
  # Covariate equations (Methods): continuous covariates use the power form
  # Y(cov) = Y * (cov / central value)^theta.  The weight central value is
  # 75 kg (stated).  The paper does not print the central values of albumin
  # or EASI ('median or another selected level of covariate'), so:
  #   * SCORE_EASI: 37, the rounded baseline EASI mean of the dupilumab arms
  #     of the same trial (37.1 for 200 mg q2w, 37.4 for 300 mg q4w; Siegfried
  #     2023 Table 1).  A trial mean, not the PK-set median.
  #   * ALB: 44 g/L, a rounded standard (neither the article nor the trial
  #     reports print a baseline albumin summary).  The same value is used by the
  #     adult sibling Kovalenko_2020_dupilumab_covariate.R, so the two
  #     remain comparable.  Albumin units are g/L (Discussion: 'SD of 3.19
  #     ... g/L').
  # A power-model centring value leaves every covariate ratio unchanged; it
  # only moves the typical value at which Vc and ke are quoted.
  #
  # Transit chain: see Kovalenko_2021_dupilumab_children_base.R (ktr = 3 / MTT;
  # three MTT boxes = depot, transit1, transit2; ka depot = transit3).
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
      notes = "Power effect on central volume Vc, (WT/75)^e_wt_vc.  Methods: the weight central value was set to 75 kg in the primary analyses for comparability across age populations.",
      source_name = "weight"
    ),
    ALB = list(
      description = "Serum albumin (baseline)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on central volume Vc, (ALB/44)^e_alb_vc.  The central value is not printed in Kovalenko 2021; 44 g/L is a rounded standard (same as the adult sibling Kovalenko_2020_dupilumab_covariate.R).  Units g/L per the Discussion (albumin SD 3.19 g/L in children).",
      source_name = "albumin"
    ),
    SCORE_EASI = list(
      description = "Eczema Area and Severity Index (baseline)",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the linear elimination rate ke, (SCORE_EASI/37)^e_score_easi_kel.  The central value is not printed in Kovalenko 2021; 37 is the rounded baseline EASI mean of the dupilumab arms of LIBERTY AD PEDS (Siegfried 2023 Am J Clin Dermatol Table 1: 37.1 and 37.4).",
      source_name = "EASI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 239L,
    n_studies = 1L,
    age_range = "6 to <12 years",
    age_mean = "8.5 (SD 1.7) years",
    weight_mean = "31.6 (SD 10.2) kg",
    sex_female_pct = 50.2,
    disease_state = "Severe atopic dermatitis, dupilumab with concomitant topical corticosteroids (LIBERTY AD PEDS, R668-AD-1652, NCT03345914).",
    dose_range = "SC dupilumab 100 mg q2w (<30 kg, n = 63) or 200 mg q2w (>=30 kg, n = 59), or 300 mg q4w (n = 119).",
    regions = "United States, Canada, Czechia, Germany, Poland, United Kingdom (Supplementary Table 1).",
    notes = "239 of 241 children on active treatment and 925 of 1173 samples in the primary analysis (Methods); 49.8% male.  BMI SD 3.35 kg/m^2 and albumin SD 3.19 g/L (Discussion).  Sparse sampling on days 1, 29, 57, 113 and end of treatment.  BLQ data handled by Beal M3.  Estimation in Monolix 2019R2 (SAEM + importance sampling)."
  )

  ini({
    # Estimated structural parameters (Table 2 / Supp. Table 3, children column)
    lvc  <- log(2.18);   label("Central volume at the covariate central values (L)")            # Table 2 children: Vc = 2.18 L (SE 0.0872)
    lkel <- log(0.0446); label("Linear elimination rate constant ke at the central EASI (1/d)")  # Table 2 children: ke = 0.0446 1/d (SE 0.00152)

    # Fixed structural parameters
    lvmax   <- fixed(log(1.64));  label("Maximum target-mediated elimination rate Vm (mg/L/d)") # Table 2 children: Vm = 1.64 (fixed; estimated in the base model per Methods)
    km      <- fixed(0.01);       label("Michaelis-Menten constant Km (mg/L)")                  # Table 2 children: Km = 0.01 (fixed)
    lkcp    <- fixed(log(0.211)); label("Central-to-peripheral rate constant kcp (1/d)")        # Table 2 children: kcp = 0.211 (fixed)
    lkpc    <- fixed(log(0.310)); label("Peripheral-to-central rate constant kpc (1/d)")        # Table 2 children: kpc = 0.310 (fixed)
    lka     <- fixed(log(0.641)); label("Absorption rate constant ka (1/d)")                    # Table 2 children: ka = 0.641 (fixed)
    lmtt    <- fixed(log(0.105)); label("Mean transit time MTT (d)")                            # Table 2 children: MTT = 0.105 (fixed)
    lfdepot <- fixed(log(0.642)); label("Subcutaneous bioavailability F (fraction)")            # Table 2 children: F = 0.642 (fixed)

    # Covariate effects (power form)
    e_wt_vc          <-  0.849; label("Power exponent of WT/75 on Vc (unitless)")                # Table 2 children: Vc ~ weight = 0.849 (SE 0.0345)
    e_alb_vc         <- -0.525; label("Power exponent of ALB/44 on Vc (unitless)")               # Table 2 children: Vc ~ albumin = -0.525 (SE 0.149)
    e_score_easi_kel <-  0.169; label("Power exponent of SCORE_EASI/37 on ke (unitless)")        # Table 2 children: ke ~ EASI = 0.169 (SE 0.0471)

    # IIV (Supp. Table 3 children): SD(ln Vc) = 0.291, SD(ln ke) = 0.417,
    # Corr(ln ke, ln Vc) = -0.883.  Variances 0.291^2, 0.417^2;
    # covariance -0.883 * 0.291 * 0.417.
    etalvc + etalkel ~ c(0.084681,
                         -0.107149, 0.173889)

    # Residual error (Supp. Table 3 children)
    propSd <- 0.131;        label("Proportional residual error (fraction)")                      # Supp. Table 3 children: sigma_prop = 13.1 CV% (SE 0.402)
    addSd  <- fixed(0.03);  label("Additive residual error (mg/L)")                              # Supp. Table 3 children: sigma_add = 0.03 mg/L (fixed)
  })
  model({
    vc   <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc * (ALB / 44)^e_alb_vc
    kel  <- exp(lkel + etalkel) * (SCORE_EASI / 37)^e_score_easi_kel
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
