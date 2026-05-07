Kovalenko_2016_dupilumab_ddmore <- function() {
  description <- paste(
    "Dupilumab population PK as encoded in DDMORE Foundation Model Repository",
    "entry DDMODEL00000273. Two-compartment model with parallel linear and",
    "Michaelis-Menten elimination from the central compartment, first-order",
    "absorption from a SC depot, and a non-standard body-weight covariate on",
    "the central volume in which weight enters log(V2) multiplicatively.",
    "Final estimates here come from the bundle's Output_simulated_*.lst (a",
    "SAEM/IMP fit on the bundle's Simulated_Dupilumab.CSV); no",
    "Output_real_*.lst is shipped, so these values do NOT match Table 2 of",
    "the publication. The publication-faithful encoding (Eq. 1, Eq. 2,",
    "Table 2 estimates) is the replicate_of counterpart at",
    "inst/modeldb/specificDrugs/Kovalenko_2016_dupilumab.R."
  )
  reference <- paste(
    "Kovalenko P, DiCioccio AT, Davis JD, Li M, Ardeleanu M, Graham NMH,",
    "Soltys R (2016). Exploratory Population PK Analysis of Dupilumab, a",
    "Fully Human Monoclonal Antibody Against IL-4Ralpha, in Atopic Dermatitis",
    "Patients and Normal Volunteers. CPT Pharmacometrics Syst Pharmacol.",
    "5(11):617-624. doi:10.1002/psp4.12136.",
    "DDMORE Foundation Model Repository: DDMODEL00000273."
  )
  vignette     <- "Kovalenko_2016_dupilumab_ddmore"
  ddmore_id    <- "DDMODEL00000273"
  replicate_of <- "inst/modeldb/specificDrugs/Kovalenko_2016_dupilumab.R"
  units        <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Bundle covariate equation (Executable_Simulated_Dupilumab.ctl $PK):",
        "C1 = exp((LWT - log(75)) * 0.75) = (WT/75)^0.75, then",
        "V2 = exp((MU_1 + ETA(1)) * C1). This makes log(V2) proportional to",
        "(WT/75)^0.75, which is mathematically distinct from the publication's",
        "Eq. 1, V2 = THETA1 * (WT/75)^THETA2. The reference weight (75 kg) is",
        "encoded in the .ctl as the constant log(75) = 4.317488."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 197L,
    n_studies      = 6L,
    age_range      = "Adults; study-wide mean age 37 years (publication)",
    age_median     = "Mean 37 years (Results, Data section of Kovalenko 2016)",
    weight_range   = "Adults; study-wide mean weight 76 kg (publication)",
    weight_median  = "Mean 76 kg (Results, Data section of Kovalenko 2016)",
    sex_female_pct = 49,
    race_ethnicity = "Not reported in the article",
    disease_state  = "Pooled healthy volunteers and patients with moderate-to-severe atopic dermatitis (publication)",
    dose_range     = "IV 1, 3, 8, 12 mg/kg single infusions; SC 75-300 mg single or weekly up to 12 doses (publication)",
    regions        = "Not specified in the article",
    notes          = paste(
      "Population demographics taken from Kovalenko 2016 (Table 1 + Results).",
      "The DDMORE bundle ships only Simulated_Dupilumab.CSV (a single-subject",
      "placeholder with WT=1 across all rows) and Output_simulated_*.lst (a",
      "SAEM/IMP fit on that simulated data); the bundle does not include the",
      "real 197-subject cohort, so the bundle's .lst final estimates are an",
      "internal regression test rather than the publication's Table 2 fit."
    )
  )

  ini({
    # Structural PK parameters - final estimates from the bundle's
    # Output_simulated_Dupilumab.lst FINAL PARAMETER ESTIMATE (SAEM block,
    # echoed by the IMP block) on the bundle's Simulated_Dupilumab.CSV.
    # These values do not match Kovalenko 2016 Table 2; the bundle's .ctl
    # parameterization differs from the publication and the simulated cohort
    # is a single-subject placeholder. See the vignette Errata.

    # V2 covariate parameter. The .ctl encodes V2 as:
    #   V2 = exp((MU_1 + ETA(1)) * (WT/75)^0.75)
    # so MU_1 (= lvc) coincides with log(V2) only at WT = 75 kg. At 75 kg
    # this gives V2 ~ 2.024 L (vs publication Table 2: 2.74 L).
    lvc      <- 0.705    ; label("V2 covariate parameter; log(V2) at 75 kg (log L)")           # .lst FINAL TH 1 (POPV2)

    # Linear elimination rate constant (1/d), log scale.  KE = exp(MU_2 + ETA(2)).
    # .lst final TH 2 = -1.40 -> typical ke = exp(-1.40) ~ 0.247 1/d (vs Table 2: 0.0459 1/d).
    lke      <- -1.40    ; label("Linear elimination rate constant ke (log 1/d)")              # .lst FINAL TH 2 (POPKE)

    # Maximum target-mediated elimination rate Vmax (mg/L/d), log scale.
    # VM = exp(MU_3 + ETA(3)).  .lst final TH 3 = -1.08 -> Vmax ~ 0.339 mg/L/d
    # (vs Table 2: 0.968 mg/L/d).
    lvmax    <- -1.08    ; label("Maximum target-mediated elimination rate Vmax (log mg/L/d)") # .lst FINAL TH 3 (POPVM)

    # First-order absorption rate constant (1/d), log scale.  KA = exp(MU_4 + ETA(4)).
    # .lst final TH 4 = -1.04 -> ka ~ 0.353 1/d (vs Table 2: 0.254 1/d).
    lka      <- -1.04    ; label("First-order absorption rate constant ka (log 1/d)")          # .lst FINAL TH 4 (POPKA)

    # Subcutaneous bioavailability.  The .ctl encodes BIO = MU_5 + ETA(5) on
    # the linear scale (no log transform); however, ETA(5) is FIXED at 0 in
    # the .ctl $OMEGA, so the linear-MU form has no IIV behaviour to preserve.
    # We store the parameter on the nlmixr2lib-canonical log scale (lfdepot)
    # for convention compliance; the back-transformed typical fraction is
    # exp(lfdepot), and the predicted concentration profile is identical to
    # the .ctl's linear-MU encoding.
    lfdepot  <- log(0.615) ; label("Subcutaneous bioavailability (log fraction)")               # .lst FINAL TH 5 (POPBIO) = 0.615 on linear scale

    # Central-to-peripheral and peripheral-to-central rate constants (1/d).
    # The .ctl encodes K23 = MU_6 + ETA(6) and K32 = MU_7 + ETA(7) on LINEAR
    # scale; the corresponding ETA(6) and ETA(7) are FIXED at 0 in the .ctl
    # $OMEGA, so there is no IIV on these.
    # .lst final TH 6 = 0.0475 (vs Table 2 k23 = 0.0652 1/d).
    # .lst final TH 7 = 0.105  (vs Table 2 k32 = 0.129 1/d).
    k23      <- 0.0475   ; label("Central-to-peripheral rate constant (1/d)")                  # .lst FINAL TH 6 (POPK23)
    k32      <- 0.105    ; label("Peripheral-to-central rate constant (1/d)")                  # .lst FINAL TH 7 (POPK32)

    # Michaelis-Menten constant (mg/L).  Fixed at 0.01 in the .ctl ($THETA
    # `(0.01 FIXED)`) per the publication, which fixed Km at 0.01 mg/L because
    # the OFV was insensitive to Km below ~0.01 mg/L.  .lst final TH 8 = 0.01.
    Km       <- fixed(0.01); label("Michaelis-Menten constant (mg/L; fixed)")                  # .lst FINAL TH 8 (POPKM)

    # Inter-individual variability.  The .ctl declares
    #   $OMEGA DIAGONAL(1) ETAV2
    #   $OMEGA BLOCK(2)    ETAKE / ETAVM (covariance structure)
    #   $OMEGA DIAGONAL(1) ETAKA
    #   $OMEGA DIAGONAL(4) FIX (ETAs 5-8)
    # The .lst SAEM final estimates of these random-effect variances are
    # essentially zero (~1e-6 to 1e-7); the SAEM run on the single-subject
    # simulated dataset has driven the IIV to a degenerate fit, so the values
    # below are NOT representative of population variability in Kovalenko 2016.
    # Use the specificDrugs/ counterpart for population-simulation work.
    etalvc                ~ 2.92e-06                                  # .lst FINAL OMEGA(1,1) (ETAV2)
    etalke + etalvmax     ~ c(2.26e-07, -1.72e-07, 1.36e-07)          # .lst FINAL OMEGA BLOCK(2): var2, cov, var3
    etalka                ~ 7.32e-08                                  # .lst FINAL OMEGA(4,4) (ETAKA)

    # Residual error.  The .ctl $ERROR computes
    #   Y = IPRE * exp(ERR(1) * SD1) + ERR(2) * SD2
    # with $SIGMA 1 FIX, so SD1 = THETA(9) is the proportional SD on the
    # log-concentration scale and SD2 = 0.03 (FIXED) is the additive SD in
    # mg/L.  Translated to nlmixr2 via the small-SD linearization
    # exp(SD1*ERR1) ~ 1 + SD1*ERR1, which is exact to within ~0.1% for SD1 of
    # the magnitude reported here.  .lst final TH 9 = 0.0377.
    propSd   <- 0.0377   ; label("Proportional residual SD on log scale")                      # .lst FINAL TH 9 (POPSD1)
    addSd    <- fixed(0.03); label("Additive residual SD (mg/L; fixed in .ctl)")               # .ctl $ERROR SD2 = 0.03
  })

  model({
    # Body-weight covariate factor as encoded in the .ctl $PK:
    #   C1 = exp((LWT - log(75)) * 0.75) = (WT/75)^0.75
    # The 0.75 exponent is hard-coded in the .ctl (not estimated) and is not
    # the same as the publication's fitted V2~weight exponent of 0.705.
    c1 <- (WT / 75)^0.75

    # Individual PK parameters
    # V2 follows the bundle's non-standard form: log(V2) is proportional to c1.
    vc   <- exp((lvc + etalvc) * c1)
    ke   <- exp(lke + etalke)
    vmax <- exp(lvmax + etalvmax)
    ka   <- exp(lka + etalka)

    # ODE system from .ctl $DES (DTYP gating dropped: route is determined by
    # the dosing compartment in the event record - SC doses to depot, IV doses
    # to central).  The .ctl's COMP=(AUC) integrator is omitted; users can
    # request AUC at solve time.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - ke * central -
                          k23 * central + k32 * peripheral1 -
                          central * vmax / (Km + central / vc)
    d/dt(peripheral1) <-  k23 * central - k32 * peripheral1

    # Subcutaneous bioavailability.  The .ctl applies F1 = BIO directly with
    # BIO = MU_5 + ETA(5) on the linear scale; ETA(5) is FIXED at 0 in the
    # .ctl, so we back-transform from the canonical log-scale parameter.
    f(depot) <- exp(lfdepot)

    # Observation: central / vc gives mg/L when dose is in mg and vc in L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
