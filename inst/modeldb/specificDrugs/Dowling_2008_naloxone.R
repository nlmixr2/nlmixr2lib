Dowling_2008_naloxone <- function() {
  description <- "Three-compartment population PK model for naloxone in six healthy male volunteers receiving 0.8 mg intravenous, 0.8 mg intramuscular, 0.8 mg intranasal, 2 mg intravenous, and 2 mg intranasal doses in an open-label crossover design (Dowling 2008). Intramuscular and intranasal absorption are modeled as first-order via separate depot compartments (Ka_im 0.65 1/h, F_im 0.36; Ka_in 1.52 1/h, F_in 0.038); intravenous doses go directly to the central compartment (F 1, structural anchor). Fat-free mass (Janmahasatian 2005 formula, called LBW2005 in the paper) is allometrically scaled on clearance with fixed exponent 0.75, and body weight is linearly scaled on central volume (exponent 1); both effects use a 70 kg reference."
  reference   <- paste(
    "Dowling J, Isbister GK, Kirkpatrick CMJ, Naidoo D, Graudins A.",
    "Population pharmacokinetics of intravenous, intramuscular, and",
    "intranasal naloxone in human volunteers.",
    "Ther Drug Monit. 2008;30(4):490-496.",
    "doi:10.1097/FTD.0b013e3181816214",
    sep = " "
  )
  vignette  <- "Dowling_2008_naloxone"
  units     <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linearly scaled on central volume V2 with fixed exponent 1 and a 70 kg reference (Results Section: V2 = TVV2 * (WT/70) * EXP(ETA2)). Per Patient Data: n=6 male volunteers, median weight 80 kg, range 75-100 kg.",
      source_name        = "WT"
    ),
    FFM = list(
      description        = "Fat-free mass computed from body weight, height, and sex via the Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065) formula.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper calls this covariate LBW2005 (lean body weight, 2005 formula) but the formula it cites is Janmahasatian's fat-free mass; nlmixr2lib canonicalises the Janmahasatian output as FFM (see inst/references/covariate-columns.md). Allometrically scaled on clearance CL with fixed exponent 0.75 and a 70 kg reference (Results Section: CL = TVCL * (LBW2005/70)^0.75 * EXP(ETA1)). Derivation for adult males: FFM = 9.27e3 * WT / (6.68e3 + 216 * BMI) with BMI = WT / height_m^2 (Janmahasatian 2005).",
      source_name        = "LBW2005"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 6L,
    n_studies      = 1L,
    age_range      = "24-45 years (median 25)",
    age_median     = "25 years",
    weight_range   = "75-100 kg (median 80)",
    weight_median  = "80 kg",
    height_range   = "1.75-1.93 m (median 1.78)",
    sex_female_pct = 0,
    race_ethnicity = NULL,
    disease_state  = "Healthy adult male volunteers with no prior or current opioid dependence, opioid analgesic use, cardiorespiratory disease, upper respiratory tract infection, or abnormal nasal anatomy.",
    dose_range     = "0.8 mg IV, 0.8 mg IM, 0.8 mg IN, 2 mg IV, and 2 mg IN across five crossover occasions per subject (minimum 2-day washout).",
    regions        = "Australia (Prince of Wales Hospital, Sydney).",
    notes          = "Open-label crossover study; 128 plasma concentrations retained above the LOQ of 1 ug/L (82 IV, 39 IM, 7 IN). Only two of six subjects had detectable concentrations above the LOQ after 2 mg IN naloxone; 0.8 mg IN produced no detectable concentrations above the LOQ (final IN dataset: two occasions). One IV 0.8 mg record was excluded for administration error. LOQ/2 imputation applied to the first sub-LOQ observation per profile (Beal M3 alternative). Race / ethnicity not reported. Sample times: 5, 10, 15, 30, 45, 60, 90, 120, 180, 240 minutes post-dose per arm. Table 1: population parameter estimates and 95% percentiles (1000 non-parametric bootstraps)."
  )

  ini({
    # Structural PK parameters (Table 1 mean estimates) -----------------------
    lka      <- log(0.65) ; label("First-order intramuscular absorption rate constant (1/h)")  # Table 1: Ka[im] = 0.65 1/h
    lka2     <- log(1.52) ; label("First-order intranasal absorption rate constant (1/h)")     # Table 1: Ka[in] = 1.52 1/h
    lcl      <- log(91)   ; label("Clearance (L/h)")                                            # Table 1: CL = 91 L/h
    lvc      <- log(2.87) ; label("Central volume of distribution V2 (L)")                     # Table 1: V2 = 2.87 L
    lvp      <- log(1.49) ; label("First peripheral volume of distribution V3 (L)")            # Table 1: V3 = 1.49 L
    lvp2     <- log(33.6) ; label("Second peripheral volume of distribution V4 (L)")           # Table 1: V4 = 33.6 L
    lq       <- log(5.66) ; label("Inter-compartmental clearance Q3 (L/h)")                    # Table 1: Q3 = 5.66 L/h
    lq2      <- log(29.8) ; label("Inter-compartmental clearance Q4 (L/h)")                    # Table 1: Q4 = 29.8 L/h
    lfdepot  <- log(0.36) ; label("Relative bioavailability of intramuscular naloxone (F_im, unitless)")  # Table 1: F_tot[im] = 0.36 (IV F fixed to 1 as structural anchor)
    lfdepot2 <- log(0.038); label("Relative bioavailability of intranasal naloxone (F_in, unitless)")    # Table 1: F_tot[in] = 0.038 (IV F fixed to 1 as structural anchor)

    # Covariate effects -------------------------------------------------------
    # Both exponents are hard-coded in the printed final-model equations
    # (Results: CL = TVCL * (LBW2005/70)^0.75; V2 = TVV2 * (WT/70)) with no
    # reported RSE / SE / CI, matching the "canonical 0.75 / 1 fixed exponent"
    # pattern (parameter-names.md, "Fixed parameters"): wrap in fixed().
    e_ffm_cl <- fixed(0.75) ; label("Allometric exponent for fat-free mass (LBW2005) on CL (unitless)")  # Results Section: CL = TVCL * (LBW2005/70)^0.75
    e_wt_vc  <- fixed(1)    ; label("Linear exponent for total body weight on V2 (unitless)")            # Results Section: V2 = TVV2 * (WT/70)^1

    # Inter-individual variability (log-normal; omega^2 = variance) -----------
    # Table 1 reports "BSV on CL (CV%)" as 0.00581 (7.6%) and "BSV on V2 (CV%)"
    # as 0.25 (50%); the numeric values are the OMEGA^2 (variance) estimates,
    # and the parenthetical percentages are approximate CV%. Verified:
    #   sqrt(exp(0.00581) - 1) = 7.63%  matches paper's 7.6%
    #   sqrt(exp(0.25)    - 1) = 53.3%  matches paper's 50% (approximation)
    # BSV was not retained on V3, V4, Q3, Q4, or Ka (Results Section).
    etalcl ~ 0.00581  # Table 1: BSV on CL, omega^2 = 0.00581 (CV ~ 7.6%)
    etalvc ~ 0.25     # Table 1: BSV on V2, omega^2 = 0.25 (CV ~ 50%)

    # Residual error ----------------------------------------------------------
    # Table 1 reports "Prop Err 0.101 (31.7%)" and "Add Err 0.001 (fixed)".
    # NONMEM SIGMA reports variances; sqrt(0.101) = 0.318 (matches 31.7% CV);
    # sqrt(0.001) = 0.0316 ug/L additive SD in ng/mL units (ng/mL == ug/L).
    propSd <- 0.318            ; label("Proportional residual error (fraction; sqrt of Table 1 SIGMA^2 = 0.101)")  # Table 1: Prop Err SIGMA^2 = 0.101 (31.7% CV)
    addSd  <- fixed(sqrt(0.001)); label("Additive residual error (ng/mL; sqrt of Table 1 SIGMA^2 = 0.001)") # Table 1: Add Err SIGMA^2 = 0.001 (fixed)
  })

  model({
    # 1. Individual parameters -----------------------------------------------
    # Covariate scaling per Results Section equations. FFM (paper's LBW2005)
    # and WT are subject-level; reference 70 kg per the printed equations.
    ka      <- exp(lka)
    ka2     <- exp(lka2)
    cl      <- exp(lcl + etalcl) * (FFM / 70)^e_ffm_cl
    vc      <- exp(lvc + etalvc) * (WT  / 70)^e_wt_vc
    vp      <- exp(lvp)
    vp2     <- exp(lvp2)
    q       <- exp(lq)
    q2      <- exp(lq2)
    fdepot  <- exp(lfdepot)
    fdepot2 <- exp(lfdepot2)

    # 2. Micro-constants -----------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system ----------------------------------------------------------
    # Compartments:
    #   depot        - intramuscular absorption depot (F_im, ka_im)
    #   depot2       - intranasal absorption depot (F_in, ka_in)
    #   central      - V2 sampling compartment (also the target for IV bolus doses)
    #   peripheral1  - V3 first peripheral compartment
    #   peripheral2  - V4 second peripheral compartment
    d/dt(depot)       <- -ka  * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <-  ka  * depot + ka2 * depot2 -
                          (kel + k12 + k13) * central +
                          k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 4. Bioavailability -----------------------------------------------------
    # IV doses (cmt = "central") bypass both depots so F = 1 implicitly.
    f(depot)  <- fdepot
    f(depot2) <- fdepot2

    # 5. Observation and error -----------------------------------------------
    # central is in mg (dose units); vc is in L; central/vc is mg/L. Multiply
    # by 1000 to express the observation in ng/mL (equivalent to ug/L, the
    # paper's assay units).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
