Kovalenko_2021_dupilumab_adolescent_covariate <- function() {
  description <- "Dupilumab primary covariate population PK model for adolescents 12 to <18 years with moderate-to-severe atopic dermatitis (Kovalenko 2021): 2-compartment model parameterised in rates with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; body weight on central volume and BMI and EASI score on the linear elimination rate."
  reference <- "Kovalenko P, Kamal MA, Davis JD, Huniti N, Xu C, Bansal A, Shumel B, DiCioccio AT. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab in Adolescents and Children >=6 to <12 Years of Age Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2021;10(11):1345-1357. doi:10.1002/cpdd.986. BMI and EASI centring values from the companion trial report: Simpson EL, et al. Efficacy and Safety of Dupilumab in Adolescents With Uncontrolled Moderate to Severe Atopic Dermatitis. JAMA Dermatol. 2020;156(1):44. doi:10.1001/jamadermatol.2019.3336, Table 1."
  vignette <- "Kovalenko_2021_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Primary COVARIATE model for adolescents >=12 to <18 years (study
  # R668-AD-1526, NCT03054428), Kovalenko 2021 Table 2 and Supplementary
  # Table 3 ('Adolescents >=12 to <18 Years of Age' column).  Same structure
  # as Kovalenko_2021_dupilumab_adolescent_base.R.
  #
  # Estimation status (Methods): Vc, ke and Vm estimated; kcp, kpc, ka, MTT,
  # F and Km fixed to the adult values.  Albumin, which is retained in the
  # children and adult covariate models, was not significant in adolescents.
  #
  # Covariate equations (Methods): continuous covariates use the power form
  # Y(cov) = Y * (cov / central value)^theta.  The weight central value is
  # 75 kg (stated).  The paper does not print the central values of BMI or
  # EASI ('median or another selected level of covariate'), so the rounded
  # baseline means of the dupilumab arms of the same trial are used
  # (Simpson 2020 Table 1; 300 mg q4w n = 84 and 200/300 mg q2w n = 82):
  #   * BMI: 24.5 kg/m^2 (arm means 24.1 and 24.9).
  #   * SCORE_EASI: 36 (arm means 35.8 and 35.3, pooled 35.55).
  # These are trial means, not the PK-set medians.  A power-model centring
  # value leaves every covariate ratio unchanged; it only moves the typical
  # value at which ke is quoted.
  #
  # Transit chain: see Kovalenko_2021_dupilumab_adolescent_base.R (ktr = 3 /
  # MTT; three MTT boxes = depot, transit1, transit2; ka depot = transit3).
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
    BMI = list(
      description = "Body mass index (baseline)",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the linear elimination rate ke, (BMI/24.5)^e_bmi_kel.  The central value is not printed in Kovalenko 2021; 24.5 kg/m^2 is the baseline BMI mean of the dupilumab arms of LIBERTY AD ADOL (Simpson 2020 JAMA Dermatol Table 1: 24.1 and 24.9).  PK-set BMI SD 6.91 kg/m^2 (Discussion).",
      source_name = "BMI"
    ),
    SCORE_EASI = list(
      description = "Eczema Area and Severity Index (baseline)",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the linear elimination rate ke, (SCORE_EASI/36)^e_score_easi_kel.  The central value is not printed in Kovalenko 2021; 36 is the rounded baseline EASI mean of the dupilumab arms of LIBERTY AD ADOL (Simpson 2020 JAMA Dermatol Table 1: 35.8 and 35.3).",
      source_name = "EASI"
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
    notes = "162 of 165 adolescents on active treatment and 827 of 1006 samples in the primary analysis (Methods); 56.7% male.  BMI SD 6.91 kg/m^2 and albumin SD 3.20 g/L (Discussion).  Sparse sampling on days 1, 29, 57, 113 and end of treatment.  BLQ data handled by Beal M3.  Estimation in Monolix 2019R2 (SAEM + importance sampling)."
  )

  ini({
    # Estimated structural parameters (Table 2 / Supp. Table 3, adolescent column)
    lvc   <- log(2.47);   label("Central volume at 75 kg (L)")                                   # Table 2 adolescents: Vc = 2.47 L (SE 0.0501)
    lkel  <- log(0.0520); label("Linear elimination rate constant ke at the central BMI and EASI (1/d)") # Table 2 adolescents: ke = 0.0520 1/d (SE 0.00188)
    lvmax <- log(1.43);   label("Maximum target-mediated elimination rate Vm (mg/L/d)")          # Table 2 adolescents: Vm = 1.43 mg/L/d (SE 0.0379)

    # Fixed structural parameters
    km      <- fixed(0.01);       label("Michaelis-Menten constant Km (mg/L)")                  # Table 2 adolescents: Km = 0.01 (fixed)
    lkcp    <- fixed(log(0.211)); label("Central-to-peripheral rate constant kcp (1/d)")        # Table 2 adolescents: kcp = 0.211 (fixed)
    lkpc    <- fixed(log(0.310)); label("Peripheral-to-central rate constant kpc (1/d)")        # Table 2 adolescents: kpc = 0.310 (fixed)
    lka     <- fixed(log(0.306)); label("Absorption rate constant ka (1/d)")                    # Table 2 adolescents: ka = 0.306 (fixed)
    lmtt    <- fixed(log(0.105)); label("Mean transit time MTT (d)")                            # Table 2 adolescents: MTT = 0.105 (fixed)
    lfdepot <- fixed(log(0.642)); label("Subcutaneous bioavailability F (fraction)")            # Table 2 adolescents: F = 0.642 (fixed)

    # Covariate effects (power form)
    e_wt_vc          <- 0.755; label("Power exponent of WT/75 on Vc (unitless)")                 # Table 2 adolescents: Vc ~ weight = 0.755 (SE 0.0517)
    e_bmi_kel        <- 0.357; label("Power exponent of BMI/24.5 on ke (unitless)")              # Table 2 adolescents: ke ~ BMI = 0.357 (SE 0.116)
    e_score_easi_kel <- 0.356; label("Power exponent of SCORE_EASI/36 on ke (unitless)")         # Table 2 adolescents: ke ~ EASI = 0.356 (SE 0.0523)

    # IIV (Supp. Table 3 adolescents): SD(ln Vc) = 0.140, SD(ln ke) = 0.304,
    # Corr(ln ke, ln Vc) = -0.529.  Variances 0.140^2, 0.304^2;
    # covariance -0.529 * 0.140 * 0.304.
    etalvc + etalkel ~ c(0.0196,
                         -0.0225142, 0.092416)

    # Residual error (Supp. Table 3 adolescents)
    propSd <- 0.0994; label("Proportional residual error (fraction)")                            # Supp. Table 3 adolescents: sigma_prop = 9.94 CV% (SE 0.602)
    addSd  <- 2.36;   label("Additive residual error (mg/L)")                                    # Supp. Table 3 adolescents: sigma_add = 2.36 mg/L (SE 0.24)
  })
  model({
    vc   <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    kel  <- exp(lkel + etalkel) * (BMI / 24.5)^e_bmi_kel * (SCORE_EASI / 36)^e_score_easi_kel
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
