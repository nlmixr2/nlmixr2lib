Kovalenko_2016_dupilumab <- function() {
  description <- "Dupilumab exploratory population PK model (Kovalenko 2016; 2-cmt with parallel linear + Michaelis-Menten elimination)"
  reference <- "Kovalenko P, DiCioccio AT, Davis JD, et al. Exploratory Population PK Analysis of Dupilumab, a Fully Human Monoclonal Antibody Against IL-4Ralpha, in Atopic Dermatitis Patients and Normal Volunteers. CPT Pharmacometrics Syst Pharmacol. 2016;5(11):617-624. doi:10.1002/psp4.12136"
  vignette <- "Kovalenko_2016_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Final model from Table 2, column "BLQ data included".  The paper reports
  # two sets of parameter estimates: "BLQ data included" (the primary model,
  # using the Beal M3 method for BLQ observations) and "BLQ data excluded"
  # (a sensitivity analysis).  The paper states "All diagnostic plots are
  # presented for the model with BLQ values" and this column is the
  # referenced final model throughout the Results section.
  #
  # Parameterization (as published): 2-compartment model with parallel linear
  # and Michaelis-Menten elimination from the central compartment, and
  # first-order absorption from a SC depot.  The paper parameterizes the
  # linear elimination as a first-order rate constant ke (1/d) acting on the
  # central amount (ke * central) rather than as clearance CL.
  # Intercompartmental transport is parameterized as k23 (1/d, central to
  # peripheral) and k32 (1/d, peripheral to central) rather than as Q and Vp;
  # the paper notes V3 = V2 * k23 / k32.  Km was fixed at 0.01 mg/L because
  # the OFV was insensitive to changes below ~0.01 mg/L; the additive
  # residual SD was fixed at 0.03 mg/L when BLQ data were included.
  #
  # IIV: Table 2 reports omega^2 (the variance of the log-scale random
  # effect) directly, so the tabulated values are inserted verbatim on the
  # right-hand side of each `~` line.
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume; reference weight 75 kg per Eq. 1 of the paper.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 197L,
    n_studies      = 6L,
    age_range      = "Adults; study-wide mean age 37 years",
    age_median     = "Mean 37 years (Results, Data section)",
    weight_range   = "Adults; study-wide mean weight 76 kg",
    weight_median  = "Mean 76 kg (Results, Data section)",
    sex_female_pct = 49,
    race_ethnicity = "Not reported in the article",
    disease_state  = "Pooled healthy volunteers and patients with moderate-to-severe atopic dermatitis",
    dose_range     = "IV 1, 3, 8, 12 mg/kg single infusions; SC 75-300 mg single or weekly up to 12 doses",
    regions        = "Not specified in the article",
    notes          = "197 participants (96 female / 101 male) across six Phase 1 and Phase 2 studies (NCT01015027, NCT01259323, NCT01385657, NCT01484600, NCT01548404, NCT01639040); 2,518 dupilumab serum concentrations analysed with NONMEM v7.3.0 and Monolix v4.3.2 using SAEM and importance sampling.  Demographic summary appears in the Results / Data subsection (page 619); study design is in Table 1."
  )

  ini({
    # Structural PK parameters - final estimates from Table 2, "BLQ data included" column
    lvc     <- log(2.74);   label("central volume at 75 kg (L)")                     # Table 2: V2 = 2.74 L
    lke     <- log(0.0459); label("linear elimination rate constant (1/d)")           # Table 2: ke = 0.0459 1/d
    lk23    <- log(0.0652); label("central-to-peripheral rate constant (1/d)")        # Table 2: k23 = 0.0652 1/d
    lk32    <- log(0.129);  label("peripheral-to-central rate constant (1/d)")        # Table 2: k32 = 0.129 1/d
    lka     <- log(0.254);  label("first-order absorption rate constant (1/d)")       # Table 2: ka = 0.254 1/d
    lvm     <- log(0.968);  label("maximum target-mediated elimination rate (mg/L/d)")# Table 2: Vm = 0.968 mg/L/d
    Km      <- fixed(0.01); label("Michaelis-Menten constant (mg/L)")                 # Table 2: Km = 0.01 (fixed)
    lfdepot <- log(0.607);  label("subcutaneous bioavailability (fraction)")          # Table 2: F = 0.607

    # Covariate effect - allometric power of weight on central volume (Eq. 1: V2 = theta1 * (WT/75)^theta2)
    e_wt_vc <- 0.705; label("weight exponent on central volume (unitless)")           # Table 2: V2 ~ weight = 0.705

    # Inter-individual variability - Table 2 reports omega^2 directly
    etalvc ~ 0.0225  # Table 2: omega^2(V2) = 0.0225
    etalke ~ 0.131   # Table 2: omega^2(ke) = 0.131
    etalka ~ 0.251   # Table 2: omega^2(ka) = 0.251
    etalvm ~ 0.0428  # Table 2: omega^2(Vm) = 0.0428

    # Residual error - Table 2: proportional CV% = 24.2 (SD = 0.242 on the linear scale);
    # additive SD fixed at 0.03 mg/L with BLQ data included.
    CcpropSd <- 0.242;       label("proportional residual error SD (fraction)") # Table 2: sigma^2 proportional (CV%) = 24.2
    CcaddSd  <- fixed(0.03); label("additive residual error SD (mg/L)")         # Table 2: sigma^2 additive = 0.03 (fixed)
  })
  model({
    # Individual PK parameters; weight-adjusted central volume per Eq. 1
    vc  <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    ke  <- exp(lke + etalke)
    k23 <- exp(lk23)
    k32 <- exp(lk32)
    ka  <- exp(lka + etalka)
    Vm  <- exp(lvm + etalvm)

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - ke * central - k23 * central + k32 * peripheral1 -
                          central * (Vm / (Km + central / vc))
    d/dt(peripheral1) <-                              k23 * central - k32 * peripheral1

    # Subcutaneous bioavailability; IV doses bypass the depot via the event record
    f(depot) <- exp(lfdepot)

    # Observation: dosing amounts in mg / central volume in L give mg/L directly
    Cc <- central / vc
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
