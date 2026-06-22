Bosch_2024_glp1ra_bodyweight <- function() {
  description <- paste(
    "QSP. GLP-1R agonist body composition model (Bosch 2024)",
    "extending the Hall 2009 three-compartment energy-balance model",
    "with (1) an inverse-Bateman lifestyle-change effect on energy",
    "intake, (2) a body-weight-dependent activity effect on physical",
    "activity energy expenditure for studies that included weight",
    "management and intensive behavioural treatment, and (3) a",
    "GLP-1R agonist drug effect driven by the in-vitro EC50-",
    "normalised free drug concentration, with a time-dependent",
    "tolerance term that shifts the in-vivo EC50 upward. Liraglutide",
    "and semaglutide PK are encoded inline as fixed one-compartment",
    "first-order absorption models (parameters from the Bosch 2024",
    "supplement S10 reproducing FDA clinical pharmacology review",
    "(liraglutide, 17 Dec 2018) and Carlsson Petri et al. 2018",
    "(semaglutide); both PK paths share the body composition",
    "machinery, and the total normalised free concentration drives",
    "the GLP-1R effect so a user simulating a single drug doses to",
    "that drug's depot only. Body weight (kg) and percent change",
    "from baseline are the primary observation outputs. Initial",
    "conditions are derived from baseline body weight, BMI, age, sex",
    "and height via the Jackson body-fat regression and the Mifflin",
    "resting-metabolic-rate equation; baseline energy intake is set",
    "to maintain steady state at PAL = 1.6 (low-active-to-active).",
    "11 active ODEs (3 macronutrient stores, 2 extracellular-water",
    "states, lipolysis diet target, adaptive thermogenesis, plus 2",
    "first-order PK chains for each drug)."
  )
  reference <- paste(
    "Bosch R, Sijbrands EJG, Snelder N. Quantification of the effect",
    "of GLP-1R agonists on body weight using in vitro efficacy",
    "information: An extension of the Hall body composition model.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13:1488-1502.",
    "doi:10.1002/psp4.13183. PMID 38867373.",
    "Body composition framework: Hall KD. Predicting metabolic",
    "adaptation, body weight change, and energy intake in humans.",
    "Am J Physiol Endocrinol Metab. 2010;298(3):E449-66.",
    "doi:10.1152/ajpendo.00559.2009.",
    "Liraglutide PK source reproduced in supplement S10:",
    "FDA Clinical Pharmacology Review for liraglutide, 17 Dec 2018.",
    "Semaglutide PK source reproduced in supplement S10:",
    "Carlsson Petri KC, Ingwersen SH, Flint A, Zacho J, Overgaard RV.",
    "Semaglutide s.c. once-weekly in type 2 diabetes: a population",
    "pharmacokinetic analysis. Diabetes Ther. 2018;9:1533-1547.",
    "doi:10.1007/s13300-018-0458-5."
  )
  vignette <- "Bosch_2024_glp1ra_bodyweight"

  paper_specific_compartments <- c(
    "fat", "prot", "carb", "decw", "bwecw",
    "lipol_diet", "therm",
    "depot_lira", "central_lira",
    "depot_sema", "central_sema"
  )

  units <- list(
    time          = "d",
    dosing        = "pmol",
    concentration = "pmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight at study entry; sets the steady-state energy balance and the allometric scaling of liraglutide and semaglutide PK.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Bosch 2024 supplement S10 uses BW0 as the baseline body weight; passed unchanged here as WT. Sets the body fat / fat-free-mass / RMR initial conditions and the (WT/85)^0.774 semaglutide CL scaling and (WT/90)^0.703 / ^1.24 liraglutide CL / Vc scaling. Hall model time-varying BW state is tracked internally as BWkg.",
      source_name        = "BW0"
    ),
    HT = list(
      description        = "Baseline body height; enters the Mifflin resting metabolic rate equation as cm and the Jackson body-fat regression via BMI.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Bosch 2024 supplement S10 carries height in metres as HGHT; the model interconverts (cm = HT, m = HT/100). RMR (kcal/d) = 9.99 * BW0 + 6.25 * HT - 4.92 * AGE + s where s = +5 for males and -161 for females.",
      source_name        = "HGHT (metres in source)"
    ),
    AGE = list(
      description        = "Subject age at baseline; enters the Mifflin RMR equation and the Jackson body-fat regression.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in both RMR (-4.92 * AGE) and Jackson body-fat regression (0.14 * AGE).",
      source_name        = "AGE"
    ),
    BMI = list(
      description        = "Baseline body mass index; drives the Jackson initial body-fat fraction via log(BMI).",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "F_init / BW0 (kg/kg) = 0.0014 * AGE + 0.3731 * log(BMI) - 1.0394 (men) or 0.0014 * AGE + 0.3996 * log(BMI) - 1.0201 (women), per Jackson 2002 Br J Nutr 89(2):277-285 as cited in supplement S10.",
      source_name        = "BMI0"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male; selects female vs male coefficients of the Jackson body-fat regression and the Mifflin RMR equation.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Bosch 2024 supplement S10 uses SEX = 0 for male, SEX = 1 for female, SEX = 99 for pooled-mean studies with a FRACFEM majority test; in nlmixr2lib the user passes a single 0/1 SEXF based on cohort majority for mean-data simulations and per-subject sex for individual-level simulations.",
      source_name        = "SEX"
    ),
    LSCI = list(
      description        = "Lifestyle change intensity: amplitude of the inverse-Bateman lifestyle-change effect on energy intake (fractional peak reduction in EI from study participation, dietary restriction, or placebo behavioural intervention).",
      units              = "fraction (0..1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Bosch 2024 supplement S10 reports study-specific LSCI values estimated per study: 0 for Can 2014 placebo, 0.079 for Pi-Sunyer 2015 SCALE, 0.099 for Astrup 2009, 0.112 for le Roux 2017 SCALE, 0.224 for STEP 1, 0.117 for STEP 1 late (TIME >= 105 d), 0.149 for STEP 5 and STEP 8 placebo arms, 0.548 for STEP 3, and -0.0288 for Blundell 2017. The user supplies the appropriate LSCI for the study/arm being simulated; the model uses it directly without rescaling.",
      source_name        = "LSCI"
    ),
    WM_IBT = list(
      description        = "Weight management with intensive behavioural treatment indicator, 1 = subject is enrolled in an arm with weight management + intensive behavioural treatment (the activity effect is applied), 0 = otherwise (no activity effect).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no WM-IBT intervention)",
      notes              = "Bosch 2024 supplement S10 gates the activity effect with an IFLAG that is set to 1 only for the STEP studies where the protocol explicitly added weight management + intensive behavioural treatment in addition to the GLP-1R agonist. The activity effect describes a study-protocol-induced increase in exercise expenditure correlated with weight loss (Bosch 2024 Figure 2); under WM_IBT = 0 the model reduces to the Hall body composition + LSC + drug-effect base.",
      source_name        = "IFLAG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 14L,
    age_range      = "18 to 75 years (pooled adult)",
    weight_range   = "67 to 117 kg (study-mean baseline)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Non-diabetic obese, pre-diabetic obese, and type-2 diabetic",
      "obese adults; diet-restriction healthy adults from the Hall",
      "2009 validation studies are also incorporated for the body",
      "composition layer."
    ),
    dose_range     = paste(
      "Liraglutide 0.6-3.0 mg SC once daily; semaglutide 0.5-2.4 mg",
      "SC once weekly (or oral semaglutide where reported). The",
      "model accepts dosing in pmol to the drug-specific depot",
      "compartment (depot_lira or depot_sema); convert from mg",
      "using the molecular weights MW_lira = 3751.20 g/mol and",
      "MW_sema = 4113.58 g/mol (e.g., 3.0 mg liraglutide = 3.0e6",
      "/ 3751.20 = 800 nmol = 800000 pmol)."
    ),
    regions        = "Global (US, EU, Asia)",
    notes          = paste(
      "Mean-study-data analysis pooling 14 publications: Diaz 1992,",
      "Jebb 1993, Jebb 1996, Schrauwen 1997, Das 2017, Heilbronn",
      "2006, Guo 2018, Redman 2007, Racette 2011, Weiss 2015,",
      "Rumpler 1991, de Boer 1986, Can 2014, Pi-Sunyer 2015 SCALE,",
      "Astrup 2009, le Roux 2017 SCALE, Blundell 2017, Hjerpsted",
      "2017, Sorli 2017 SUSTAIN 1, Wadden 2021 STEP 3, Garvey 2022",
      "STEP 5, Rubino 2022 STEP 8, Pratley 2018 SUSTAIN 7. STEP 1",
      "(Wilding 2021) was reserved as external validation. The",
      "primary fittable endpoint is body weight (kg) and percent",
      "change from baseline; the model also captures fat mass and",
      "fat-free mass internally for diagnostic plots. The Kred rate",
      "of return-to-baseline of the lifestyle effect was estimated",
      "as 0.00195 / d for non-STEP drug studies, 0.00541 / d",
      "(0.0379 / week) for STEP 1 / 5 / 8, and 0.00924 / d (0.0647",
      "/ week) for STEP 3; the model exposes a single Kred ini()",
      "parameter set to the non-STEP value (0.00195 / d) and lets a",
      "user override per simulation, or estimate per-study Kred in",
      "a downstream re-fit."
    )
  )

  ini({
    # ===== GLP-1R agonist drug effect (Bosch 2024 Table 3 Model D) =====
    ltvglp50  <- log(48.1)
    label("Log of TVGLP50 = in-vivo / in-vitro EC50 fold (unitless); time-invariant component of the in-vivo EC50")  # Bosch 2024 Table 3 Model D (48.1, 25.2% RSE)
    lt50      <- log(1439)
    label("Log of T50 = time at which the in-vivo EC50 has shifted halfway to its steady-state plateau (days)")     # Bosch 2024 Table 3 Model D (1439 d, 40.3% RSE)
    lssec50   <- log(3798)
    label("Log of SSEC50 = steady-state additive shift in in-vivo EC50 (unitless, in the same normalised units as TVGLP50)") # Bosch 2024 Table 3 Model D (3798, 32.2% RSE)
    lfulira   <- log(0.00243)
    label("Log of free fraction of liraglutide (unitless)")                                                         # Bosch 2024 Table 3 Model D (0.00243, 4.19% RSE); ~99.76 percent plasma protein binding

    # ===== Activity effect on physical activity energy expenditure =====
    lemaxact  <- log(0.00193)
    label("Log of EMAX of the activity effect on the exercise part of PAE (kcal/d/kg)")                             # Bosch 2024 Table 3 Model D (0.00193 kcal/d/kg)
    ewc50act  <- -22.9
    label("Weight change at half-maximal activity effect (kg, signed - negative because weight loss is the active direction)") # Bosch 2024 Table 3 Model D (-22.9 kg)

    # ===== Lifestyle-change effect =====
    lkdiet    <- fixed(log(10))
    label("Log of Kdiet = inverse-Bateman onset rate constant for the LSC effect (1/d, FIXED to a large value so the LSC effect onset is rapid)") # Bosch 2024 supplement S10 (THETA5 = 10, FIX)
    lkred     <- log(0.00195)
    label("Log of Kred = inverse-Bateman reduction rate constant for the LSC effect (1/d) - default non-STEP-study value") # Bosch 2024 Table 3 Model D (Kred = 0.00195 / d for non-STEP studies; STEP 1/5/8 = 0.00541 / d, STEP 3 = 0.00924 / d - user may override per study)

    # ===== Compound in-vitro potency and free fraction (Bosch 2024 Table 2, FIXED) =====
    lec50lira <- fixed(log(1.2))
    label("Log of in-vitro liraglutide EC50 (pM, FIXED) - GLP-1R activation assay")                                 # Bosch 2024 Table 2 (1.2 pM, in absence of serum protein)
    lec50sema <- fixed(log(0.9))
    label("Log of in-vitro semaglutide EC50 (pM, FIXED) - GLP-1R activation assay")                                 # Bosch 2024 Table 2 (0.9 pM, in absence of serum protein)
    lfusema   <- fixed(log(0.0025))
    label("Log of semaglutide free fraction (unitless, FIXED)")                                                     # Bosch 2024 Table 2 (0.25%)

    # ===== Liraglutide PK (Bosch 2024 supplement S10; values from FDA Clinical Pharmacology Review, 17 Dec 2018) =====
    lkalira     <- fixed(log(0.0608 * 24))
    label("Log of liraglutide first-order absorption rate constant (1/d, FIXED)")                                    # supplement S10 (KAlira = 0.0608 * 24 = 1.4592 / d)
    lcllira_ref <- fixed(log(1.11 * 24))
    label("Log of liraglutide clearance at reference covariates (L/d, FIXED; reference 90 kg male)")                 # supplement S10 (CLlira = 1.11 * 24 = 26.64 L/d at WT = 90 kg, male)
    lvclira_ref <- fixed(log(0.16))
    label("Log of liraglutide central volume of distribution at reference covariates (L, FIXED; reference 90 kg male)") # supplement S10 (VClira = 0.16 L at WT = 90 kg, male)
    e_wt_cl_lira   <- fixed(0.703)
    label("Allometric exponent on WT for liraglutide CL (unitless, FIXED)")                                          # supplement S10 ((BW0 / 90)^0.703)
    e_wt_vc_lira   <- fixed(1.24)
    label("Allometric exponent on WT for liraglutide Vc (unitless, FIXED)")                                          # supplement S10 ((BW0 / 90)^1.24)
    e_male_cl_lira <- fixed(log(1.32))
    label("Log multiplier on liraglutide CL for males (unitless, FIXED; applied as multiplier 1.32 if SEXF = 0)")    # supplement S10 (1.32^MALE)
    e_male_vc_lira <- fixed(log(1.4))
    label("Log multiplier on liraglutide Vc for males (unitless, FIXED; applied as multiplier 1.4 if SEXF = 0)")     # supplement S10 (1.4^MALE)

    # ===== Semaglutide PK (Bosch 2024 supplement S10; values from Carlsson Petri et al. 2018) =====
    lkasema     <- fixed(log(0.0286 * 24))
    label("Log of semaglutide first-order absorption rate constant (1/d, FIXED)")                                    # supplement S10 (KAsema = 0.0286 * 24 = 0.6864 / d)
    lvcsema_ref <- fixed(log(12.2))
    label("Log of semaglutide central volume of distribution at reference covariates (L, FIXED; not WT-scaled in supplement)") # supplement S10 (VCsema = 12.2 L)
    lclsema_ref <- fixed(log(0.0478 * 24))
    label("Log of semaglutide clearance at reference covariates (L/d, FIXED; reference 85 kg)")                       # supplement S10 (CLsema = 0.0478 * 24 = 1.1472 L/d at WT = 85 kg)
    e_wt_cl_sema <- fixed(0.774)
    label("Allometric exponent on WT for semaglutide CL (unitless, FIXED)")                                          # supplement S10 ((BW0 / 85)^0.774)

    # ===== Residual error =====
    propSd    <- 0.32
    label("Proportional residual error on body weight (BWkg) - approximate CV from supplement S10 SIGMA block")        # Bosch 2024 supplement S10 ($SIGMA @block 0.105221 0 0.209354; sqrt(0.105) approx 0.325)
  })

  model({
    # -- Drug effect ini-block back-transforms --
    tvglp50  <- exp(ltvglp50)
    t50      <- exp(lt50)
    ssec50   <- exp(lssec50)
    fulira   <- exp(lfulira)
    fusema   <- exp(lfusema)
    ec50lira <- exp(lec50lira)
    ec50sema <- exp(lec50sema)
    emaxact  <- exp(lemaxact)

    # -- Lifestyle change rate constants --
    kdiet <- exp(lkdiet)
    kred  <- exp(lkred)

    # -- PK back-transforms with covariate adjustments --
    # Liraglutide: MALE = 1 - SEXF; supplement multiplier 1.32^MALE on CL, 1.4^MALE on Vc.
    male_indicator <- 1 - SEXF
    kalira <- exp(lkalira)
    cllira <- exp(lcllira_ref + e_wt_cl_lira * log(WT / 90) + e_male_cl_lira * male_indicator)
    vclira <- exp(lvclira_ref + e_wt_vc_lira * log(WT / 90) + e_male_vc_lira * male_indicator)
    kellira <- cllira / vclira

    # Semaglutide: WT scaling on CL only (Vc not scaled per supplement).
    kasema <- exp(lkasema)
    vcsema <- exp(lvcsema_ref)
    clsema <- exp(lclsema_ref + e_wt_cl_sema * log(WT / 85))
    kelsema <- clsema / vcsema

    # ===== Initial-condition derivations (Bosch 2024 supplement S10) =====
    # Jackson body-fat regression - male equation gives the male reference.
    # Body fat (g) = (BW0 / 100) * (0.14 * AGE + 37.31 * log(BMI) - 103.94) * 1000.
    # SEXF = 1 uses the female coefficients (39.96, 102.01).
    f_frac_male   <- (0.14 * AGE + 37.31 * log(BMI) - 103.94) / 100
    f_frac_female <- (0.14 * AGE + 39.96 * log(BMI) - 102.01) / 100
    f_frac        <- f_frac_male * (1 - SEXF) + f_frac_female * SEXF
    finit_g       <- WT * f_frac * 1000
    bwinit_g      <- WT * 1000
    ffminit_g     <- bwinit_g - finit_g

    # Mifflin RMR equation - kcal/d
    ht_cm <- HT
    rmr_male   <- 9.99 * WT + 6.25 * ht_cm - 4.92 * AGE + 5
    rmr_female <- 9.99 * WT + 6.25 * ht_cm - 4.92 * AGE - 161
    rmr_init   <- rmr_male * (1 - SEXF) + rmr_female * SEXF
    pal_init   <- 1.6
    ei_init    <- rmr_init * pal_init

    # Macronutrient caloric densities (kcal/g)
    rho_f <- 9.44
    rho_p <- 4.7
    rho_c <- 4.18
    rho_k <- 4.45

    # Digestibility (Southgate & Durnin 1970)
    digest_c <- 0.95
    digest_f <- 0.96
    digest_p <- 0.90

    # Standard meal composition (Hall 2009)
    f_carb <- 0.5
    f_prot <- 0.13
    f_fat  <- 0.37

    # Baseline caloric intake split (kcal/d)
    ci0 <- f_carb * ei_init
    pi0 <- f_prot * ei_init
    fi0 <- f_fat  * ei_init

    # Molecular and stoichiometric masses
    g_mass        <- 180
    tg_mass       <- 860
    glycerol_mass <- 92
    ffa_mass      <- (rho_f * tg_mass - rho_c * glycerol_mass) / (3 * rho_f)
    n_mass        <- 14
    k_mass        <- 104
    aa_mass       <- 110

    glycerol_exog <- glycerol_mass / (rho_f * tg_mass)
    atp_kcal      <- 19
    atp_tg_synth  <- 8
    atp_synth_p   <- 5
    atp_proteol_p <- 1
    atp_g_synth   <- 2
    atp_nexcr     <- 4

    # Stoichiometric energy efficiencies
    dnleff <- 0.37 * rho_f / rho_c
    gngeff <- 0.8
    keteff <- 4.5 * rho_k * k_mass / (rho_f * ffa_mass)
    eta_f  <- atp_tg_synth * atp_kcal / tg_mass
    eta_p  <- atp_synth_p  * atp_kcal / aa_mass
    eps_p  <- atp_proteol_p * atp_kcal / aa_mass
    eta_g  <- atp_g_synth  * atp_kcal / g_mass

    # Body composition (Minnesota Experiment / Wang AJCN 2003)
    pfrac_cm   <- 0.25
    icw_frac_cm <- 0.7
    hg         <- 2.7
    hp         <- 1.6

    bw_keys   <- 67533
    f_keys    <- 9050
    g_keys    <- 500
    ffm_keys  <- bw_keys - f_keys
    bm_keys   <- 0.04 * bw_keys
    ecw_keys  <- 0.7 * 0.235 * bw_keys
    ecp_keys  <- 0.732 * bm_keys + 0.01087 * ecw_keys
    p_keys    <- pfrac_cm * (ffm_keys - bm_keys - ecp_keys - ecw_keys)

    bm        <- 0.04 * WT * 1000
    ecw_b     <- 0.7 * 0.235 * WT * 1000
    ecp       <- 0.732 * bm + 0.01087 * ecw_b
    ginit     <- 500
    cm_init   <- ffminit_g - bm - ecp - ecw_b
    pinit     <- pfrac_cm * cm_init
    kicw      <- (icw_frac_cm - pfrac_cm * hp) * cm_init - ginit * hg
    ics       <- (1 - pfrac_cm - icw_frac_cm) * cm_init - ginit
    icw_b     <- kicw + pinit * hp + ginit * hg

    # SMR slopes / values
    smr_liver  <- 0.2
    smr_skmusc <- 0.013
    smr_brain  <- 0.24
    smr_kid    <- 0.44
    smr_hrt    <- 0.44
    smr_res    <- 0.012
    slope_liver  <- 0.01736
    slope_skmusc <- 0.5934
    slope_kid    <- 0.003786
    slope_hrt    <- 0.00288
    slope_res    <- 0.3747
    smr_ffmb <- smr_liver  * slope_liver  + smr_skmusc * slope_skmusc +
                smr_kid    * slope_kid    + smr_hrt    * slope_hrt    +
                smr_res    * slope_res
    smr_f <- 0.0045
    mbrain <- 1400

    # TEF (thermic effect of food) fractions
    tef_f <- 0.025
    tef_c <- 0.075
    tef_p <- 0.25

    # Baseline metabolizable intakes (kcal/d, post-digestibility, in kcal/d as the ODE substrate units)
    cin_b <- digest_c * ci0
    pin_b <- digest_p * pi0
    fin_b <- digest_f * fi0
    ei_b  <- cin_b + pin_b + fin_b
    tef_init <- tef_f * fin_b + tef_c * cin_b + tef_p * pin_b

    # Hill / DNL initialisation
    hill_dnl <- 4
    kdnl     <- 2
    dnl_init <- cin_b * (ginit / g_keys)^hill_dnl /
                (kdnl^hill_dnl + (ginit / g_keys)^hill_dnl)

    # Lipolysis baseline
    lipol_b      <- 0.16
    s1s          <- 2.0 / 3.0
    lipol_init   <- (finit_g / f_keys)^s1s
    glycerolprod_init <- rho_c * lipol_b * lipol_init * glycerol_mass
    gngf_b_endog <- lipol_b * glycerol_mass * rho_c
    gngf_b_exog  <- glycerol_exog * fin_b * rho_c
    gngf_init    <- gngf_b_endog * lipol_init + gngf_b_exog

    # Ketogenesis
    maxkg_p     <- 80
    maxkg_np    <- 40
    maxkg_npg   <- 20
    avekg_p     <- 10
    betak_p     <- log(maxkg_p / maxkg_np)
    betak_g     <- log(maxkg_np / maxkg_npg)
    k_klip      <- maxkg_npg / avekg_p - 1
    ketogen_pct_init <- maxkg_p * exp(-betak_g * ginit / g_keys) *
      exp(-betak_p) * lipol_init / (k_klip + lipol_init)
    ketogen_init <- tg_mass * lipol_b * lipol_init * ketogen_pct_init / 100
    ketox_init   <- rho_k * ketogen_init

    # Proteolysis / GNG baselines
    proteol_b <- 2.73
    proteol_init <- proteol_b * (pinit / p_keys)
    gngp_b <- 300
    gngp_init <- gngp_b * (proteol_init / proteol_b)

    # Glycogenolysis baseline
    glycogenol_b <- 1
    glycogenol_init <- glycogenol_b * (ginit / g_keys)
    gmin <- 10

    # Substrate-competition fitted parameters (Hall 2010 Table 1)
    activ_vs_rest <- 0.52
    wci_over_wg   <- 0.93
    wpi           <- 1.1
    spi_neg       <- 1.7
    sci           <- 0.85
    spi_plus      <- 3.8
    tau_sig       <- 1.1
    tau_lip       <- 1 / log(2)
    pa_lip        <- 0.4
    sldiet        <- 2
    kl_diet       <- 4
    sg            <- 1
    lipol_60hr    <- 3.1
    lipol_min     <- 0.9
    lipol_max     <- lipol_60hr / (1 - exp(-2.5 / tau_lip))
    k_lip         <- log((lipol_max - lipol_min) / (1 - lipol_min))

    therm_const_uf <- 0.74
    therm_const_of <- 0.02
    tau_therm      <- 7

    # ECW dynamics
    na_conc       <- 140 * 23 / 1000
    na_zero_cin   <- 4000
    na_ecw        <- 3
    tau_ecw_bw    <- 1000
    ecw_inc       <- 0.16
    na_b          <- 4000

    # GNG carbohydrate / protein cross-coupling
    gng_ci <- 4 * 6.25 * rho_p / gngp_b
    gng_pi <- (0.56 - 0.2 * gng_ci) / 1.5
    glyc_gng_effect <- 0
    diet_pto <- 0

    # Baseline exercise term (kcal/d/kg)
    exerc_b <- 0

    # ===== Time-varying state algebra =====
    ecw       <- ecw_b + decw

    # Free drug normalised concentrations (sum across both drugs;
    # only the dosed drug has non-zero state in single-drug simulations).
    c_lira    <- central_lira / vclira
    c_sema    <- central_sema / vcsema
    cdrugf_norm <- c_lira * fulira / ec50lira + c_sema * fusema / ec50sema

    # GLP-1R drug effect with time-dependent tolerance (Bosch 2024 Equations 7a-b)
    glp50     <- tvglp50 + (ssec50 * t) / (t50 + t)
    glpeff    <- 1 - cdrugf_norm / (glp50 + cdrugf_norm)

    # Lifestyle change effect on EI (Bosch 2024 Equation 5).
    # The inverse-Bateman equals 1 at t = 0 (exp(0) - exp(0) = 0), so
    # no special-case for t = 0 is required.
    lsceff    <- 1 - (LSCI * kdiet * (exp(-kdiet * t) - exp(-kred * t))) / (kred - kdiet)

    # Modified intakes
    cin <- digest_c * (ci0 * lsceff * glpeff)
    pin <- digest_p * (pi0 * lsceff * glpeff)
    fin <- digest_f * (fi0 * lsceff * glpeff)

    dcin <- cin - cin_b
    dpin <- pin - pin_b
    dfin <- fin - fin_b
    ei   <- cin + pin + fin

    # Sodium balance
    na_imbal <- na_b - na_b - na_ecw * (ecw - ecw_b) - na_zero_cin * (1 - cin / cin_b)

    # Body composition assembly
    cm   <- kicw + prot * (1 + hp) + carb * (1 + hg) + ics
    ffm  <- bm + ecp + ecw + cm
    bw_g <- fat + ffm
    bwkg <- bw_g / 1000

    # Organ masses
    mliver  <- 588   + slope_liver  * (bm + ecp + ecw_b + cm)
    mskmusc <- -6301 + slope_skmusc * (bm + ecp + ecw_b + cm)
    mkid    <- 134   + slope_kid    * (bm + ecp + ecw_b + cm)
    mhrt    <- 1.5 * (7.06 + slope_hrt * (bm + ecp + ecw_b + cm))
    mres    <- 4818  + slope_res    * (bm + ecp + ecw_b + cm)

    # Baseline activity coefficient (Hall 2010); the exact closed-form
    # baseline-correction Ec below uses activ_b, so it is defined first.
    activ_b <- -exerc_b / 1000 +
      (ei_b - tef_init - rmr_init) / bwinit_g

    # Baseline-correction energy term Ec (Bosch 2024 supplement S10).
    # Ec is the constant offset that closes the books of the energy
    # balance at baseline so that TEE = EI_b at t = 0.
    ec <- ei_b - (tef_init + (activ_b + exerc_b / 1000) * bwinit_g +
                  smr_ffmb * (ffminit_g - mbrain) +
                  smr_brain * mbrain +
                  smr_f * finit_g +
                  (1 - dnleff) * dnl_init +
                  rho_k * (1 - keteff) * ketogen_init +
                  (1 - gngeff) * (gngp_init + gngf_init) +
                  (eta_p + eps_p) * proteol_init * aa_mass +
                  eta_f * lipol_b * lipol_init * tg_mass +
                  eta_g * glycogenol_init * g_mass +
                  ((atp_nexcr * atp_kcal / n_mass) / (rho_p * 6.25)) * pin_b)

    # Activity effect (Bosch 2024 Equations 6a-c) - gated by WM_IBT covariate
    wtch_kg <- bwkg - WT
    acteff_raw <- emaxact * (wtch_kg / (ewc50act + wtch_kg))
    acteff <- WM_IBT * acteff_raw
    exerc  <- exerc_b + acteff

    # Adaptive thermogenesis split (algebraic, from current Therm state)
    dactiv_eff  <- activ_vs_rest * therm
    dsmrl_eff   <- (1 - activ_vs_rest) * therm
    smr_ffm     <- smr_ffmb * (1 + dsmrl_eff)
    activ_tsb   <- activ_b
    activ_eff   <- activ_tsb * (1 + dactiv_eff)

    # Thermogenesis driver depends on feed status: under-feeding uses
    # therm_const_uf, over-feeding uses therm_const_of.
    therm_const <- therm_const_of + (therm_const_uf - therm_const_of) * (ei < ei_b)

    # Lipolysis (g/d) - Hall 2010
    lipol_act    <- pa_lip * ((activ_tsb + exerc) / (activ_b + exerc_b / 1000) - 1)
    lipol_state  <- (fat / f_keys)^s1s * (lipol_diet + lipol_act)
    mfat         <- fat / f_keys - 1
    mfat2        <- max(0.0, mfat)
    lipol_diet_target <- 1 + ((lipol_max - lipol_min) * exp(-k_lip * cin / cin_b) +
                              lipol_min - 1) *
                              kl_diet^sldiet / (mfat2^sldiet + kl_diet^sldiet)
    lipol_rate   <- lipol_b * lipol_state * tg_mass

    # Ketogenesis (g/d) - downstream of lipolysis
    ketogen_pct <- maxkg_p * exp(-betak_g * carb / g_keys) *
      exp(-betak_p * pin / pin_b) * lipol_state / (k_klip + lipol_state)
    ketogen   <- lipol_rate * ketogen_pct / 100
    kgmax     <- 400
    kspill    <- 70
    kumax     <- 20
    kurine_raw <- ketogen * kumax / (kgmax - kspill) - kumax / (kgmax / kspill - 1)
    kurine    <- (ketogen >= kspill) * kurine_raw
    ketox_g   <- ketogen - kurine
    ketox     <- ketox_g * rho_k

    # Nutrient balance constants (Hall 2010 Eq. A.10 family)
    kcbal_num   <- cin_b - dnl_init
    kcbal_denom <- ei_b - (1 - dnleff) * dnl_init - rho_k * (1 - keteff) * ketogen_init -
      gngp_init - gngf_init + glycerolprod_init - ketox_init
    kcbal <- kcbal_num / kcbal_denom

    kfbal_num   <- 3 * ffa_mass / tg_mass * fin_b + dnleff * dnl_init -
      rho_k * (1 - keteff) * ketogen_init - ketox_init
    kfbal <- kfbal_num / kcbal_denom

    kpbal_num   <- pin_b - gngp_init
    kpbal <- kpbal_num / kcbal_denom

    wg  <- (g_keys / ginit)^sg * (kcbal / kpbal) * (pinit / p_keys + wpi) /
      (1 + (g_keys / ginit)^sg * wci_over_wg * ginit / (gmin + ginit))
    wf  <- kfbal / (1 - kfbal) * (1 + kcbal / kpbal) *
      (pinit / p_keys + wpi) * (f_keys / finit_g)^s1s
    wci <- wg * wci_over_wg

    # Proteolysis (g/d)
    proteol_state <- prot / p_keys + diet_pto * (pin - pin_b) / pin_b
    proteol_rate  <- proteol_b * proteol_state * aa_mass

    # Gluconeogenesis from protein
    tcarb <- tanh(carb / g_keys - 1)
    mtcarb <- 1 - glyc_gng_effect * tcarb
    gngp_glycogen_effect <- max(0.0, mtcarb)
    gng1 <- gngp_b * (proteol_state + gng_pi * dpin / pin_b -
                       gng_ci * dcin / cin_b) * gngp_glycogen_effect
    gngp <- max(0.0, gng1)
    gngf_endog <- gngf_b_endog * lipol_state
    gngf_exog  <- glycerol_exog * fin * rho_c
    gngf       <- gngf_endog + gngf_exog

    # Glycogenolysis (g/d)
    glycogenol_state <- max(0.0, carb / g_keys)
    glycogenol_rate  <- glycogenol_b * glycogenol_state * g_mass

    # Inactivity / protein-intake constraints (Hall 2010 Activ_max derivation)
    nbal_inactive <- -2
    delta_energy_inact <- -(activ_b + exerc_b / 1000) * bwinit_g -
      (atp_nexcr * atp_kcal / n_mass) * nbal_inactive
    delta_g_inactiv <- -delta_energy_inact / (2 * rho_c)
    inact_p <- (pin_b - gngp_init - rho_p * 6.25 * nbal_inactive) /
      (ei_b + delta_energy_inact - (1 - dnleff) * dnl_init -
       rho_k * (1 - keteff) * ketogen_init - gngp_init - gngf_init +
       glycerolprod_init - ketox_init)
    wg_inactiv <- wg * (1 + delta_g_inactiv / ginit)
    pi_inact <- (p_keys / pinit) * ((inact_p / (1 - inact_p)) *
      (wci + wg_inactiv + (1 - pa_lip) * wf * (finit_g / f_keys)^s1s) - wpi)
    activ_max <- max(1.0, pi_inact)
    activ_min <- 0

    # Substrate competition activity factor
    k_act <- log((activ_max - activ_min) / (1 - activ_min))
    activity_factor <- (activ_max - activ_min) *
      exp(-k_act * (activ_tsb + exerc) / (activ_b + exerc_b / 1000)) + activ_min

    spi <- spi_plus + (spi_neg - spi_plus) * (dpin < 0)

    # psig is reserved as an additional state in the supplement but
    # never integrated (DADT(8) = 0 in active runs); held at zero here
    # so wpi2 collapses to wpi.
    psig <- 0
    wf2  <- wf * lipol_state
    fat_term <- max(0.0, wf2)
    wci2 <- wci * (1 + sci * dcin / cin_b)
    mwci2 <- max(0.0, wci2)
    carb_term <- wg * glycogenol_state^sg + mwci2 * carb / (gmin + carb)
    wpi2 <- wpi * (1 + psig)
    prot_term <- max(0.0, wpi2) + proteol_state * activity_factor

    z_total <- carb_term + prot_term + fat_term
    carbfrac <- carb_term / z_total
    protfrac <- prot_term / z_total
    fatfrac  <- fat_term  / z_total

    # DNL during dynamics
    dnl <- cin * (carb / g_keys)^hill_dnl /
      (kdnl^hill_dnl + (carb / g_keys)^hill_dnl)

    # TEF
    tef <- tef_f * fin + tef_c * cin + tef_p * pin

    # Energy expenditure backbone (Hall 2010 Eq. A.18-A.20). Ec closes
    # the baseline-energy balance at t = 0 (TEE = EI_b).
    ebar_back <- ec + (activ_eff + exerc) * bw_g +
      smr_ffm * (ffm - mbrain) -
      smr_ffm * (ecw - ecw_b) -
      smr_ffm * (carb - g_keys) * (1 + hg) +
      smr_brain * mbrain +
      smr_f * fat +
      (1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen +
      (1 - gngeff) * (gngp + gngf) +
      (eta_p + eps_p) * proteol_rate +
      eta_f * lipol_rate + eta_g * glycogenol_rate

    nuke <- (atp_nexcr * atp_kcal / n_mass) / (rho_p * 6.25)

    # Stoichiometric correction factor in TEE denominator
    denom <- 1 + (eta_p / rho_p - nuke) * protfrac +
      (eta_f / rho_f) * fatfrac +
      (eta_g / rho_c) * carbfrac

    # GlycerolProd algebra (Bosch 2024 supplement S10 closed-form).
    # Step 1: build EbarB - the energy expenditure backbone *without*
    # GlycerolProd in the nuke term (GlycerolProd is not yet known here).
    nonfat_init <- (1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen
    ebarB <- ebar_back + tef +
      nuke * (gngp * (1 - protfrac) - protfrac *
                (nonfat_init + gngf + ketox))
    # Step 2: NumeratorA, the TEE-like quantity used to derive GlycerolProd.
    num_a <- ebarB +
      (eta_g / rho_c) * (cin - dnl + carbfrac *
                          (nonfat_init + gngp + gngf + ketox)) +
      (eta_p / rho_p) * (pin - gngp + protfrac *
                          (nonfat_init + gngp + gngf + ketox)) +
      (eta_f / rho_f) * (3 * ffa_mass / tg_mass * fin + dnleff * dnl -
                          rho_k * (1 - keteff) * ketogen - rho_k * kurine -
                          ketox + fatfrac *
                          (nonfat_init + gngp + gngf + ketox))
    g_corr <- rho_c * glycerol_mass / (rho_f * tg_mass)
    glycerolprod <- ((lipol_rate / tg_mass) * rho_c * glycerol_mass +
                     g_corr * (3 * ffa_mass / tg_mass * fin + dnleff * dnl -
                               rho_k * (1 - keteff) * ketogen -
                               rho_k * kurine - ketox -
                               fatfrac * (num_a / denom -
                                          (1 - dnleff) * dnl -
                                          rho_k * (1 - keteff) * ketogen -
                                          gngp - gngf - ketox))) /
                    (1 + g_corr * fatfrac / denom)

    # Step 3: NumeratorB - the full numerator with GlycerolProd included.
    num_b <- ebar_back +
      tef +
      nuke * (gngp * (1 - protfrac) - protfrac *
                ((1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen + gngf -
                 glycerolprod + ketox)) +
      (eta_g / rho_c) * (cin - dnl + carbfrac *
                          ((1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen +
                           gngp + gngf - glycerolprod + ketox)) +
      (eta_p / rho_p) * (pin - gngp + protfrac *
                          ((1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen +
                           gngp + gngf - glycerolprod + ketox)) +
      (eta_f / rho_f) * (3 * ffa_mass / tg_mass * fin + dnleff * dnl -
                          rho_k * (1 - keteff) * ketogen - rho_k * kurine -
                          ketox + fatfrac *
                          ((1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen +
                           gngp + gngf - glycerolprod + ketox))
    tee <- num_b / denom

    # Macronutrient ODEs (Hall 2010 Eq. A.21-A.23)
    nonfat_term <- (1 - dnleff) * dnl + rho_k * (1 - keteff) * ketogen
    d_dt_fat_g  <- (3 * ffa_mass / tg_mass * fin + dnleff * dnl -
                    rho_k * (1 - keteff) * ketogen - rho_k * kurine - ketox -
                    fatfrac * (tee - nonfat_term - gngp - gngf + glycerolprod - ketox)) /
                    rho_f
    d_dt_prot_g <- (pin - gngp -
                    protfrac * (tee - nonfat_term - gngp - gngf + glycerolprod - ketox)) /
                    rho_p
    d_dt_carb_g <- (cin - dnl -
                    carbfrac * (tee - nonfat_term - gngp - gngf + glycerolprod - ketox)) /
                    rho_c

    # ===== Differential equations =====
    d/dt(fat)  <- d_dt_fat_g
    d/dt(prot) <- d_dt_prot_g
    d/dt(carb) <- d_dt_carb_g
    d/dt(decw) <- (1 / na_conc) * na_imbal + bwecw
    d/dt(bwecw) <- (ecw_inc * (bw_g - bwinit_g) - bwecw) / tau_ecw_bw
    d/dt(lipol_diet) <- (lipol_diet_target - lipol_diet) / tau_lip
    d/dt(therm) <- (therm_const * (ei - ei_b) / ei_b - therm) / tau_therm

    # Drug PK ODEs - first-order absorption, one-compartment disposition,
    # one chain per drug; user doses to depot_lira OR depot_sema.
    d/dt(depot_lira)   <- -kalira * depot_lira
    d/dt(central_lira) <-  kalira * depot_lira - kellira * central_lira
    d/dt(depot_sema)   <- -kasema * depot_sema
    d/dt(central_sema) <-  kasema * depot_sema - kelsema * central_sema

    # Initial conditions
    fat(0)         <- finit_g
    prot(0)        <- pinit
    carb(0)        <- ginit
    decw(0)        <- 0
    bwecw(0)       <- 0
    lipol_diet(0)  <- 1
    therm(0)       <- 0

    # ===== Observations =====
    Cc        <- bwkg
    wtch_pct  <- 100 * (bwkg - WT) / WT
    Cc ~ prop(propSd)
  })
}
