Kovalenko_2021_dupilumab_children_base <- function() {
  description <- "Dupilumab primary base population PK model for children 6 to <12 years with severe atopic dermatitis (Kovalenko 2021): 2-compartment model parameterised in rates with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; body weight (reference 75 kg) is the only covariate, on central volume."
  reference <- "Kovalenko P, Kamal MA, Davis JD, Huniti N, Xu C, Bansal A, Shumel B, DiCioccio AT. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab in Adolescents and Children >=6 to <12 Years of Age Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2021;10(11):1345-1357. doi:10.1002/cpdd.986"
  vignette <- "Kovalenko_2021_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Primary BASE model for children >=6 to <12 years (study R668-AD-1652,
  # NCT03345914), Kovalenko 2021 Table 1 and Supplementary Table 2 ('Children
  # >=6 to <12 Years of Age' column).  The structure is the adult model of
  # Kovalenko 2020 (doi:10.1002/cpdd.780) applied without structural change
  # (Methods, 'Population PK Analysis'; Figure 1).
  #
  # Estimation status (Methods and Table 1 footnote):
  #   * Vc, ke and Vm estimated on the phase 3 data.  Table 1 prints Vm as
  #     '1.64 (fixed)' for children, but its footnote gives the SE of Vm
  #     'in the pediatric base model, where it was estimated rather than
  #     fixed as 0.0511 mg/L/d', and the Methods state 'Vm was estimated
  #     using the base model and fixed in the covariate model'.  Vm is
  #     therefore estimated (not fixed()) here and fixed in the covariate
  #     file.
  #   * ka = 0.641 1/d was estimated on the semisparse phase 2a
  #     R668-AD-1412 data and then fixed (Methods); it is the only structural
  #     value that differs from the adult fixed set.
  #   * kcp, kpc, MTT, F and Km fixed to the adult values.
  #
  # Transit chain: Figure 1 shows the SC dose entering a chain of three
  # transit compartments spanned by the mean transit time (MTT), which feed
  # the absorption depot drained at ka.  The three MTT boxes are depot,
  # transit1 and transit2 here, and the ka depot is transit3; ktr = 3 / MTT
  # so that MTT is the mean time spent in the chain.  The same NN = 3,
  # KTR = NN / MTT convention is printed in the deposited control stream of
  # the later dupilumab analysis by the same group (Nguyen 2026,
  # doi:10.1002/cpt.70233; see Nguyen_2026_dupilumab.R).
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
    n_subjects = 239L,
    n_studies = 1L,
    age_range = "6 to <12 years",
    age_mean = "8.5 (SD 1.7) years",
    weight_mean = "31.6 (SD 10.2) kg",
    sex_female_pct = 50.2,
    disease_state = "Severe atopic dermatitis, dupilumab with concomitant topical corticosteroids (LIBERTY AD PEDS, R668-AD-1652, NCT03345914).",
    dose_range = "SC dupilumab 100 mg q2w (<30 kg, n = 63) or 200 mg q2w (>=30 kg, n = 59), or 300 mg q4w (n = 119).",
    regions = "United States, Canada, Czechia, Germany, Poland, United Kingdom (Supplementary Table 1).",
    notes = "239 of 241 children on active treatment and 925 of 1173 samples in the primary analysis (Methods); 49.8% male.  Sparse sampling on days 1, 29, 57, 113 and end of treatment.  BLQ data handled by Beal M3.  Estimation in Monolix 2019R2 (SAEM + importance sampling)."
  )

  ini({
    # Estimated structural parameters (Table 1 / Supp. Table 2, children column)
    lvc   <- log(2.22);   label("Central volume at 75 kg (L)")                                   # Table 1 children: Vc = 2.22 L (SE 0.0945)
    lkel  <- log(0.0444); label("Linear elimination rate constant ke (1/d)")                     # Table 1 children: ke = 0.0444 1/d (SE 0.00155)
    lvmax <- log(1.64);   label("Maximum target-mediated elimination rate Vm (mg/L/d)")          # Table 1 children: Vm = 1.64 mg/L/d; footnote: estimated in the base model, SE 0.0511

    # Fixed structural parameters
    km      <- fixed(0.01);       label("Michaelis-Menten constant Km (mg/L)")                  # Table 1 children: Km = 0.01 (fixed)
    lkcp    <- fixed(log(0.211)); label("Central-to-peripheral rate constant kcp (1/d)")        # Table 1 children: kcp = 0.211 (fixed)
    lkpc    <- fixed(log(0.310)); label("Peripheral-to-central rate constant kpc (1/d)")        # Table 1 children: kpc = 0.310 (fixed)
    lka     <- fixed(log(0.641)); label("Absorption rate constant ka (1/d)")                    # Table 1 children: ka = 0.641 (fixed; estimated on phase 2a R668-AD-1412 data per Methods)
    lmtt    <- fixed(log(0.105)); label("Mean transit time MTT (d)")                            # Table 1 children: MTT = 0.105 (fixed)
    lfdepot <- fixed(log(0.642)); label("Subcutaneous bioavailability F (fraction)")            # Table 1 children: F = 0.642 (fixed)

    # Covariate effect
    e_wt_vc <- 0.864; label("Power exponent of WT/75 on Vc (unitless)")                          # Table 1 children: Vc ~ weight = 0.864 (SE 0.0371)

    # IIV (Supp. Table 2 children): SD(ln Vc) = 0.305, SD(ln ke) = 0.409,
    # Corr(ln ke, ln Vc) = -0.871.  Variances 0.305^2, 0.409^2;
    # covariance -0.871 * 0.305 * 0.409.
    etalvc + etalkel ~ c(0.093025,
                         -0.108653, 0.167281)

    # Residual error (Supp. Table 2 children)
    propSd <- 0.132;        label("Proportional residual error (fraction)")                      # Supp. Table 2 children: sigma_prop = 13.2 CV% (SE 0.402)
    addSd  <- fixed(0.03);  label("Additive residual error (mg/L)")                              # Supp. Table 2 children: sigma_add = 0.03 mg/L (fixed)
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
