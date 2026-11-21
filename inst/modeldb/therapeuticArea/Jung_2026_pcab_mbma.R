Jung_2026_pcab_mbma <- function() {
  description <- paste(
    "MBMA. Class-level semi-mechanistic PK/PD model for potassium-competitive",
    "acid blockers (PCABs: tegoprazan, YH4808, fexuprazan, vonoprazan) fitted",
    "to digitised arm-mean plasma concentrations and intragastric pH from 13",
    "published healthy-volunteer studies (Jung 2026). One-compartment PK with",
    "first-order absorption whose rate is scaled by a Gaussian function of the",
    "current intragastric pH; an indirect-response pH state whose production",
    "rate kin is modulated by a 24-h cosine circadian term and by a tanh-bounded",
    "stimulus summing a power function of the central drug AMOUNT and the amount",
    "of food in a first-order gastric-emptying food compartment. Each model",
    "individual is a study ARM (digitised mean profile), so the per-ID random",
    "effects are between-arm variability; drug identity is a super-level",
    "(NONMEM $LEVEL) random effect on ka, CL, V, the drug-stimulus scaling",
    "factor and the pH of maximal absorption, nested on DIDN. Only the",
    "across-drug centre is published: the per-drug estimates of ka, V and the",
    "drug-stimulus scaling factor are not reported, so the model simulates a",
    "typical or randomly drawn PCAB, not any one named drug. The amount entered",
    "per meal into the food compartment is not reported, so the food-effect",
    "scaling factor has no stated unit and meal dosing is left to the user."
  )
  reference <- paste(
    "Jung W, Lee J, Jeon H, Sung T, Yun H, Lee S, Chae JW.",
    "A Mechanism-Based Multi-Level Population PK/PD Model for",
    "Potassium-Competitive Acid Blockers.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70181.",
    "doi:10.1002/psp4.70181"
  )
  vignette <- "Jung_2026_pcab"
  units <- list(
    time = "h",
    dosing = "mg (drug into depot); meal amount into food in user-chosen units (not reported in the source)",
    concentration = "ng/mL (Cc); gastric_ph is in pH units"
  )

  paper_specific_compartments <- "food"

  compartmentData <- list(
    depot = list(analyte = "PCAB", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "PCAB", units = "mg", specimen = "plasma", verified = TRUE),
    gastric_ph = list(analyte = "intragastric pH", units = "pH units", specimen = "not applicable", verified = TRUE),
    food = list(
      analyte = "food in the stomach",
      units = "(meal amount, units not reported)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    DIDN = list(
      description = "Integer drug-identity index used as the super-level nesting column for the inter-drug random effects on ka, CL, V, the drug-stimulus scaling factor and the pH of maximal absorption.",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- DIDN is a nesting level, not a covariate with a reference category. It indexes ETA(9)-ETA(13), declared in the Data S2 control stream as $LEVEL DRUG=(9:13[4]).",
      notes = "Jung 2026 coded DRUG = 1 tegoprazan, 2 YH4808, 3 fexuprazan, 4 vonoprazan (Results, Data Collection). The drug-level empirical Bayes estimates are reported only for SCALE (Discussion), so a DIDN value does NOT select a named drug's parameters: each distinct DIDN value receives one random draw from the inter-drug distribution, shared by every arm (ID) carrying that value. Consumed by the ini() nesting syntax, so it is not referenced in model(). Simulation caveat: an event table must carry at least two distinct DIDN values, and omega must be passed explicitly (omega = mod$omega) because a nested omega is a list of matrices keyed by level.",
      source_name = "DRUG"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected where reported (Jung 2026 Methods) but not tested or retained: 'due to incomplete reporting across studies, these variables were not incorporated as statistical covariates'."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Collected where reported (Jung 2026 Methods) but not tested or retained, for the same reason as AGE."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Collected where reported (Jung 2026 Methods; control-stream column SEX) but not tested or retained, for the same reason as AGE."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 507L,
    n_studies = 13L,
    n_arms = 50L,
    age_range = "not reported (healthy adult volunteers)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported; source studies conducted in Korea, Japan, the UK and China",
    disease_state = "Healthy volunteers.",
    dose_range = "Single and multiple oral doses: tegoprazan 50-400 mg, YH4808 30-600 mg, fexuprazan 10-320 mg, vonoprazan 1-120 mg (Jung 2026 Table 1).",
    regions = "Korea, Japan, UK, China",
    n_observations = "PK and PD observations per arm in Jung 2026 Table 1.",
    notes = "Digitised arm-mean plasma concentration and intragastric pH profiles from 13 published studies (Jung 2026 Table 1: Chung 2022 YH4808; Sunwoo 2018 fexuprazan; Sakurai 2015 and Jenkins 2015 vonoprazan; Han 2019, Hwang 2019, Ghim 2021, He 2021, Yoon 2021, Jeon 2021, Han 2021, Sunwoo 2020, Yang 2022 tegoprazan). Each of the 50 study arms is one model individual; n_subjects = 507 is the sum of the per-arm subject counts in Table 1 (9 YH4808 arms with 24 subjects, 10 fexuprazan arms with 80, 14 vonoprazan arms with 176, 17 tegoprazan arms with 227). Drug identity is a super-level grouping column."
  )

  ini({
    # Final estimates: Jung 2026 Table 2 'Estimated model parameters'.
    # Structure: Data S2 control stream. Every structural theta is
    # exponentiated in the code ($THETAI THETA(3:15) = LOG(THETAI(3:15))),
    # so each is carried on the log scale here.

    # ---- Physiology of intragastric pH ----
    lkin <- log(0.435)
    label("Zero-order production rate of intragastric pH, kin (1/h)") # Table 2 No. 3 'Kin (1/h)' = 0.435 (RSE 20.9%)
    lrbase <- log(1.090)
    label("Baseline intragastric pH, Base (pH units)") # Table 2 No. 4 'Base' = 1.090 (RSE 3.2%)
    lemax <- log(5.570)
    label("Maximal stimulus-driven increase in steady-state pH, pHMax (pH units)") # Table 2 No. 5 'pHMax' = 5.570 (RSE 3.2%)
    lkst <- log(1.820)
    label("First-order gastric emptying rate constant of food, Kg (1/h)") # Table 2 No. 6 'K g (1/h)' = 1.820 (RSE 20.1%)
    lamp <- log(0.704)
    label("Amplitude of the 24-h circadian pH rhythm, AMP (pH units)") # Table 2 No. 7 'AMP' = 0.704 (RSE 13.8%)

    # ---- pH feedback on absorption (Gaussian) ----
    lph_ka_opt <- log(4.080)
    label("Intragastric pH at which absorption is fastest, SCALE (pH units)") # Table 2 No. 8 'SCALE' = 4.080 (RSE 6.7%)
    lph_ka_width <- log(2.660)
    label("Width of the Gaussian pH-absorption function, SHAPE (pH units)") # Table 2 No. 9 'SHAPE' = 2.660 (RSE 8.9%)

    # ---- PK ----
    lka <- log(0.156)
    label("Absorption rate constant at the optimal pH, Ka (1/h)") # Table 2 No. 10 'K a (1/h)' = 0.156 (RSE 7.4%)
    lcl <- log(98.700)
    label("Apparent clearance, CL (L/h)") # Table 2 No. 11 'CL (1/h)' = 98.700 (RSE 5.9%); the printed unit is a typo for L/h (control stream: 'clearance (L/h)')
    lvc <- log(56.600)
    label("Apparent volume of distribution, Vd (L)") # Table 2 No. 12 'V d (L)' = 56.600 (RSE 17.8%)

    # ---- Stimulus on pH production ----
    lstim_food <- log(0.043)
    label("Food-stimulus scaling factor, Sf (per unit meal amount)") # Table 2 No. 13 'S f' = 0.043 (RSE 51.1%)
    lpow_food <- fixed(log(1))
    label("Power exponent on the food amount, Slopef (unitless)") # Table 2 No. 14 'Slopef' = 1 (fixed); the control stream writes SFOOD*A(4) with no exponent
    lstim_central <- log(0.693)
    label("Drug-stimulus scaling factor on the central amount, Sd (1/mg)") # Table 2 No. 15 'S d (1/mg)' = 0.693 (RSE 18.3%)
    lpow_central <- log(0.572)
    label("Power exponent on the central drug amount, Sloped (unitless)") # Table 2 No. 16 'Sloped' = 0.572 (RSE 15.2%)

    # ---- Between-arm variability (ID level; Table 2 'IIV, CV%') ----
    # omega^2 = log(CV^2 + 1); e.g. CL IDV 161.14% gives 1.280, matching the
    # Data S2 $OMEGA initial 1.28431. The control stream carries KIN, BASE and
    # AMP as an OMEGA BLOCK(3) whose final covariances are not reported, so the
    # three are diagonal here.
    eta_arm_lkin ~ 0.67595 # Table 2 No. 3 Kin IIV 98.28% CV (RSE 42.9%, shrinkage 15%)
    eta_arm_lrbase ~ 0.0059883 # Table 2 No. 4 Base IIV 7.75% CV (RSE 150.5%, shrinkage 15%)
    eta_arm_lamp ~ 0.30793 # Table 2 No. 7 AMP IIV 60.05% CV (RSE 40.5%, shrinkage 26%)
    eta_arm_lph_ka_opt ~ 0.081985 # Table 2 No. 8 SCALE IIV 29.23% CV (RSE 31.8%, shrinkage 17%)
    eta_arm_lka ~ 0.040969 # Table 2 No. 10 Ka IIV 20.45% CV (RSE 38.3%, shrinkage 35%)
    eta_arm_lcl ~ 0.083988 # Table 2 No. 11 CL IIV 29.60% CV (RSE 30.3%, shrinkage 20%)
    eta_arm_lvc ~ 0.076991 # Table 2 No. 12 Vd IIV 28.29% CV (RSE 42.3%, shrinkage 56%)
    eta_arm_lstim_central ~ 0.26395 # Table 2 No. 15 Sd IIV 54.96% CV (RSE 39.4%, shrinkage 23%)

    # ---- Inter-drug variability (super level; Table 2 'IDV, CV%') ----
    # Data S2: $LEVEL DRUG=(9:13[4]), LEVCENTER=1. Table 2 prints 31.56% for
    # both the Ka and the Sd IDV, with different RSEs; both carried as printed.
    etalph_ka_opt_drug ~ 0.011984 | DIDN # Table 2 No. 8 SCALE IDV 10.98% CV (RSE 58.2%)
    etalka_drug ~ 0.094950 | DIDN # Table 2 No. 10 Ka IDV 31.56% CV (RSE 124.3%)
    etalcl_drug ~ 1.2800 | DIDN # Table 2 No. 11 CL IDV 161.14% CV (RSE 9.5%)
    etalvc_drug ~ 0.58696 | DIDN # Table 2 No. 12 Vd IDV 89.36% CV (RSE 50.6%)
    etalstim_central_drug ~ 0.094950 | DIDN # Table 2 No. 15 Sd IDV 31.56% CV (RSE 66.7%)

    # ---- Residual error ----
    expSd <- 0.545
    label("Log-scale residual SD on plasma concentration (log-additive)") # Table 2 No. 1 'Prop err (PK)' = 0.545 (RSE 4.5%); Data S2 Y = LOG(CONC) + THETA(1)*EPS(1)
    addSd_gastric_ph <- 0.880
    label("Additive residual SD on intragastric pH (pH units)") # Table 2 No. 2 'Add err (PD)' = 0.880 (RSE 4.2%); Data S2 W = THETA(2)/PH on LOG(PH), the log-scale form of an additive pH error
  })

  model({
    # Physiology (Data S2 $PK items 1-4)
    kin <- exp(lkin + eta_arm_lkin)
    rbase <- exp(lrbase + eta_arm_lrbase)
    kout <- kin / rbase
    emax <- exp(lemax)
    kst <- exp(lkst)
    amp <- exp(lamp + eta_arm_lamp)

    # Gaussian pH feedback on absorption (Data S2 $PK item 5)
    ph_ka_opt <- exp(lph_ka_opt + etalph_ka_opt_drug + eta_arm_lph_ka_opt)
    ph_ka_width <- exp(lph_ka_width)

    # Drug-specific parameters: drug-level and arm-level etas (Data S2 $PK item 6)
    ka <- exp(lka + etalka_drug + eta_arm_lka)
    cl <- exp(lcl + etalcl_drug + eta_arm_lcl)
    vc <- exp(lvc + etalvc_drug + eta_arm_lvc)
    stim_central <- exp(lstim_central + etalstim_central_drug + eta_arm_lstim_central)
    pow_central <- exp(lpow_central)
    stim_food <- exp(lstim_food)
    pow_food <- exp(lpow_food)
    kel <- cl / vc

    # Circadian term on kin, normalised by the baseline. t is the clock hour
    # of the day (t = 0 is the circadian pH peak), per the Methods ('Time is
    # the hour of the day').
    circ <- amp * cos(2 * pi * t / 24) / rbase

    # Stimulus: power of the central drug AMOUNT (mg), not concentration,
    # plus the food amount; Data S2 guards the power with IF(A(2)>0).
    if (central > 0) {
      stim_drug <- stim_central * central^pow_central
    } else {
      stim_drug <- 0
    }
    stim_meal <- 0
    if (food > 0) {
      stim_meal <- stim_food * food^pow_food
    }
    stim <- stim_drug + stim_meal
    # Data S2: PHS = PHMAX*(2/(1+EXP(-2*SCON)) - 1)/BASE, i.e. tanh(SCON).
    # The Methods print PHS = pHmax x 1/(1+Scon), which would lower pH as the
    # stimulus rises; the executable code is followed.
    phs <- emax * (2 / (1 + exp(-2 * stim)) - 1) / rbase

    feed <- exp(-((gastric_ph - ph_ka_opt)^2) / ph_ka_width^2)

    d/dt(depot) <- -ka * feed * depot
    d/dt(central) <- ka * feed * depot - kel * central
    d/dt(gastric_ph) <- (1 + circ + phs) * kin - kout * gastric_ph
    d/dt(food) <- -kst * food

    # Data S2: A_0(3) = BASE + AMP
    gastric_ph(0) <- rbase + amp

    # Central amount in mg over V in L is mg/L; x 1000 gives ng/mL (Data S2
    # S2 = V2/1000).
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
    gastric_ph ~ add(addSd_gastric_ph)
  })
}
