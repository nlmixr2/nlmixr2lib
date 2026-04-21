Kovalenko_2020_dupilumab <- function() {
  description <- "Dupilumab PK model (Kovalenko 2020)"
  reference <- "Kovalenko P, Davis JD, Li M, et al. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2020;9(6):756-767. doi:10.1002/cpdd.780"
  vignette <- "Kovalenko_2020_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Model 1 from Table 1 and Supplementary Table S2 of Kovalenko 2020.
  #
  # Parameterization (as published): 2-compartment model with parallel linear +
  # Michaelis-Menten elimination from central, and a 3-transit-compartment
  # absorption model off a SC depot.  The paper parameterizes the linear
  # elimination as a first-order rate constant ke (1/d) acting on the central
  # amount (ke * central), rather than as clearance CL.  Intercompartmental
  # transport is parameterized as kcp (1/d) and kpc (1/d) directly (with
  # Mpc = kcp/kpc), rather than as Q and Vp.  Km and F were fixed in Model 1 to
  # values carried forward from Kovalenko 2016 (doi:10.1002/psp4.12136); the
  # file reproduces those fixings.
  #
  # IIV: the paper explicitly defines omega as the standard deviation (SD) of
  # the between-subject random effect (see Methods).  The nlmixr2 `~` syntax on
  # the RHS of an eta line stores the VARIANCE (omega^2), so the published SDs
  # are squared in the ini() block below.
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume; reference weight 75 kg (assumed from prior Kovalenko 2016 publication; not explicitly stated in Kovalenko 2020).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 2041L,
    n_studies      = 16L,
    age_range      = "Adults and adolescents; detailed age breakdown not reported in the main text",
    age_median     = "Not reported in the main text",
    weight_range   = "Not reported in the main text",
    weight_median  = "Not reported in the main text (Kovalenko 2016 used 75 kg as the allometric reference weight)",
    sex_female_pct = "Not reported in the main text",
    race_ethnicity = "Race was a tested covariate; detailed breakdown not reported in the main text (primarily White, with Black, Asian, and Other categories represented across Phase 3 AD trials).",
    disease_state  = "Moderate-to-severe atopic dermatitis in adults and adolescents, plus healthy volunteers (202 HV / 1913 AD).",
    dose_range     = "IV and SC dosing pooled across 16 studies; the approved adult maintenance regimen is 300 mg Q2W SC (with 600 mg SC loading dose per the label).",
    regions        = "Multi-regional Phase 1-3 programme; specific regional breakdown not reported in the main text.",
    notes          = "Total pooled cohort: N = 2115 (2041 on active treatment, 18,243 of 20,809 samples analysed). Source studies listed in Table S1 include R668-AS-0907, TDU12265, PKM14161, R668-AD-1117 (Model 1), several Phase 1/2 studies for Model 2, and Phase 3 studies R668-AD-1334 (SOLO 1), R668-AD-1416 (SOLO 2), and R668-AD-1224 (CHRONOS). Main-text tables report structural parameters only; detailed demographics are not reproduced in the published article."
  )

  ini({
    lvc <- log(2.48); label("central volume (L)")
    lke <- log(0.0534); label("elimination rate (1/d)")
    lkcp <- log(0.213); label("central-to-peripheral rate (1/d)")
    Mpc <- 0.686; label("ratio of kcp and kpc (kpc is peripheral to central rate with units of 1/d)")
    lka <- log(0.256); label("absorption rate (1/d)")
    lmtt <- log(0.105); label("mean transit time (d)")
    lvm <- log(1.07); label("maximum target-mediated rate of elimination (mg/L/d)")
    Km <- fixed(0.01); label("Michaelis-Menten constant (mg/L)")
    lfdepot <- log(0.643); label("Bioavailability (fraction)")
    e_wt_vc <- 0.711; label("Exponent of weight on central volume (unitless)")

    # Kovalenko 2020 explicitly defines omega as the SD (standard deviation) of
    # between-subject variability ("omega (omega, standard deviation [SD] of
    # between-subject variability)").  nlmixr2's `etalxxx ~ value` syntax stores
    # the VARIANCE (omega^2) of the random effect, so the paper's SD values are
    # squared here.  The point estimates (0.192, 0.285, 0.474, 0.236, 0.525) are
    # reported as SDs in Supplementary Table 2.
    etalvc  ~ 0.192^2  # Supp. Table S2: omega_Vc  (SD) = 0.192
    etalke  ~ 0.285^2  # Supp. Table S2: omega_ke  (SD) = 0.285
    etalka  ~ 0.474^2  # Supp. Table S2: omega_ka  (SD) = 0.474
    etalvm  ~ 0.236^2  # Supp. Table S2: omega_Vm  (SD) = 0.236
    etalmtt ~ 0.525^2  # Supp. Table S2: omega_MTT (SD) = 0.525; applied on log(MTT) here to prevent negative MTT draws (a reparameterization of the Supp. Table 2 additive-on-MTT formulation)

    CcpropSd <- 0.15; label("Proportional residual error (fraction)")
    CcaddSd <- fixed(0.03); label("Additive residual error (mg/L)")
  })
  model({
    # Weight normalization to 75 kg is assumed based on prior publication.  It
    # is not specified in the current publication:
    # Kovalenko P, DiCioccio AT, Davis JD, et al. Exploratory Population PK
    # Analysis of Dupilumab, a Fully Human Monoclonal Antibody Against
    # IL-4Ralpha, in Atopic Dermatitis Patients and Normal Volunteers. CPT
    # Pharmacometrics Syst Pharmacol. 2016;5(11):617-624. doi:10.1002/psp4.12136
    vc <- exp(lvc + etalvc)*(WT/75)^e_wt_vc
    ke <- exp(lke + etalke)
    kcp <- exp(lkcp)
    ka <- exp(lka + etalka)
    MTT <- exp(lmtt + etalmtt)
    Vm <- exp(lvm + etalvm)

    # Derived parameters
    kpc <- kcp/Mpc
    ktr <- (3 + 1)/MTT

    d/dt(depot) <- -ktr*depot
    d/dt(transit1) <- ktr*(depot - transit1)
    d/dt(transit2) <- ktr*(transit1 - transit2)
    d/dt(transit3) <- ktr*transit2 - ka*transit3
    # Linear and Michaelis-Menten clearance
    d/dt(central) <-                 ka*transit3 - ke*central - kcp*central + kpc*periph - central*(Vm/(Km + central/vc))
    d/dt(periph) <-                                             kcp*central - kpc*periph

    f(depot) <- exp(lfdepot)
    # No unit conversion is required to change mg/L (dosing amount/central
    # volume unit) to mg/L (measurement unit)
    Cc <- central/vc
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
