Back_2018_fenofibrate <- function() {
  description <- paste(
    "Mechanism-based oral absorption / disposition model for fenofibrate",
    "(parent) and fenofibric acid (active form, measured analyte) in",
    "healthy Korean adults under fasted, standard-meal, and high-fat-meal",
    "conditions. Three drug compartments (stomach -> duodenum -> central)",
    "coupled to a 2-compartment calorie sub-model (stomach -> duodenum)",
    "via a bile-acid-driven coupling: the combined fenofibrate-metabolism /",
    "fenofibric-acid-absorption rate constant km&a is multiplied by",
    "(1 + Ebile * calories_in_duodenum), and a time-varying gastric",
    "emptying rate constant kg is multiplied by (1 + Efood) for the first",
    "6.94 h after a meal. Meal-type-specific shifts on Vc/F encode the",
    "additional bioavailability change between fasted, standard, and",
    "high-fat meals."
  )
  reference <- paste(
    "Back H, Song B, Pradhan S, Chae J, Han N, Kang W, Chang MJ, Zheng J,",
    "Kwon K, Karlsson MO, Yun H. A mechanism-based pharmacokinetic model of",
    "fenofibrate for explaining increased drug absorption after food",
    "consumption. BMC Pharmacology and Toxicology. 2018;19:5.",
    "doi:10.1186/s40360-018-0194-5"
  )
  vignette <- "Back_2018_fenofibrate"
  paper_specific_compartments <- c("stomach_food", "duodenum_food")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Per-dose-record covariate. Set FED = 1 for the standard-meal and",
        "high-fat-meal arms of the Back 2018 three-way crossover study and",
        "FED = 0 for the fasted arm. Drives (a) the additive Vc/F shift",
        "via the e_fed_vc_std and e_fed_vc_hf coefficients (in combination",
        "with FED_HIGHFAT to discriminate standard vs high-fat meal) and",
        "(b) the time-varying gastric-emptying boost on kg during the",
        "first 6.94 h after dosing."
      ),
      source_name        = "Food / Fed indicator"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat meal at dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-high-fat meal: fasted or standard meal)",
      notes              = paste(
        "Per-dose-record covariate refining FED. Set FED_HIGHFAT = 1 for",
        "the high-fat-meal arm (1280 kcal, 35.5% fat) and FED_HIGHFAT = 0",
        "for the fasted and standard-meal (686.3 kcal, 19.9% fat) arms.",
        "Selects between the two meal-type Vc/F coefficients in Back 2018",
        "Table 3 (E_Vc1 = -0.394 for standard meal, E_Vc2 = -0.461 for",
        "high-fat meal). Standard-meal indicator inside model() is",
        "FED * (1 - FED_HIGHFAT); high-fat-meal indicator is FED_HIGHFAT."
      ),
      source_name        = "Meal type (Standard meal / High-fat meal)"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued occasion / period indicator for inter-occasion",
        "variability (IOV). Values 1 = fasted, 2 = standard meal,",
        "3 = high-fat meal in the Back 2018 three-way crossover design."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into binary indicators oc1, oc2, oc3",
        "and multiplexed across per-occasion etas on lvc and lkel. The",
        "IOV variance for Vc/F (50.9% CV) and kel (44.9% CV) is",
        "single-valued per parameter -- occasion-2 and occasion-3 etas",
        "are fixed equal to occasion-1 (NONMEM $OMEGA BLOCK(1) SAME).",
        "Set OCC consistently across all observations within a study",
        "period; OCC indexes the period, not the order in which periods",
        "were assigned per subject."
      ),
      source_name        = "OCC / Period"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 24L,
    n_studies       = 1L,
    age_range       = "mean 23 years",
    weight_range    = "mean 68.75 kg",
    height_range    = "mean 173.29 cm",
    sex_female_pct  = 45.8,
    race_ethnicity  = "Korean (24/24)",
    disease_state   = "Healthy volunteers",
    dose_range      = "250 mg fenofibrate (SR) capsule, single oral dose, three-way crossover (fasted, standard meal, high-fat meal)",
    regions         = "South Korea (Chungnam National University)",
    n_observations  = NA_integer_,
    notes           = paste(
      "Three-way crossover food-effect study conducted April-November 2002.",
      "Each subject received a single 250 mg fenofibrate SR capsule with",
      "240 mL water 10 min after consuming no food, a standard breakfast,",
      "or a high-fat breakfast (see Back 2018 Table 1 for meal",
      "composition). Plasma samples for fenofibric acid (the active",
      "metabolite measured as a proxy for fenofibrate due to rapid",
      "esterase-mediated conversion) collected at predose and 1, 2, 3, 4,",
      "5, 6, 8, 10, 12, 24, 48, and 72 h post-dose. NONMEM 7.3 FOCE-I,",
      "1000-replicate bootstrap (Back 2018 Methods 'Pharmacokinetic",
      "modeling' and 'Model evaluation'). Physiologic GI volumes fixed",
      "from Physiology of the Gastrointestinal Tract (Elsevier, 2006",
      "ref [19]): V_stomach = 49 mL fasted / 1 L fed, V_duodenum = 45 mL;",
      "these are reported in Table 3 but do not enter the rate equations",
      "of the amount-based ODE system implemented here."
    )
  )

  ini({
    # ============================================================
    # Structural parameters. Estimated values from Back 2018 Table 3
    # ('Estimated parameters from the final MBPK model') final-model
    # column. Time units are 1/h; volume is L.
    # ============================================================
    lkg <- log(0.0412); label("Fasted gastric emptying rate constant kg (1/h)")                        # Table 3 final: kg = 0.0412 1/h (RSE 8.5%)
    lkma <- log(0.198); label("Combined fenofibrate-metabolism / fenofibric-acid-absorption rate km&a (1/h)") # Table 3 final: km&a = 0.198 1/h (RSE 28.9%)
    lkel <- log(0.27);  label("Fenofibric-acid elimination rate constant kel (1/h)")                   # Table 3 final: kel = 0.27 1/h (RSE 13.6%)
    lvc <- log(12.9);   label("Apparent central volume of distribution of fenofibric acid Vc/F (L)")   # Table 3 final: Vc/F = 12.9 L (RSE 27.8%)
    lkgp <- log(0.00971); label("Stomach-to-duodenum calorie transit rate constant kg' (1/h)")          # Table 3 final: kg' = 0.00971 1/h (RSE 58.9%)
    lkout <- log(0.00972); label("Duodenum calorie elimination rate constant kout (1/h)")               # Table 3 final: kout = 0.00972 1/h (RSE 38.3%)

    # ============================================================
    # Food / bile-acid effects on absorption.
    # ============================================================
    e_food_kg <- 0.617;  label("Fed-state multiplicative boost on kg during the first mtime2 hours after a meal (unitless)") # Table 3 final: E_food = 0.617 (RSE 40.4%)
    e_cal_kma <- 0.0239; label("Bile-acid-mediated per-calorie multiplicative effect on km&a (1/kcal)") # Table 3 final: E_bile = 0.0239 (RSE 22.2%)
    mtime2    <- fixed(6.94); label("Duration of the post-meal kg boost (h; held at the estimate during simulation)") # Table 3 final: MTIME2 = 6.94 h (RSE 12.4%); held fixed in simulation as a structural break-point

    # ============================================================
    # Meal-type-specific Vc/F shifts (additive coefficients in
    # 1 + E_Vc1*standard + E_Vc2*highfat parameterisation; both
    # coefficients are negative -> higher F under fed conditions).
    # ============================================================
    e_fed_vc_std <- -0.394; label("Standard-meal multiplicative offset on Vc/F (applied as (1 + e_fed_vc_std) when FED * (1 - FED_HIGHFAT) = 1)") # Table 3 final: E_Vc1 = -0.394 (RSE 40.1%)
    e_fed_vc_hf  <- -0.461; label("High-fat-meal multiplicative offset on Vc/F (applied as (1 + e_fed_vc_hf) when FED_HIGHFAT = 1)")              # Table 3 final: E_Vc2 = -0.461 (RSE 37.7%)

    # ============================================================
    # Inter-individual variability. Back 2018 Table 3 reports IIV as
    # %CV; the stored log-normal variance is omega^2 = log(1 + CV^2).
    # kg and kg' share a single IIV (Back 2018 Table 3 footnote
    # "k_g (k_g', shared IIV)" reports a single 31.7% IIV column for
    # both rate constants).
    # ============================================================
    etalkg ~ 0.09575                                              # Table 3 IIV kg (shared with kg'): 31.7% CV  -> omega^2 = log(1 + 0.317^2)
    etalkel ~ 0.55660                                             # Table 3 IIV kel:  86.3% CV -> omega^2 = log(1 + 0.863^2)
    etalvc  ~ 0.62317                                             # Table 3 IIV Vc/F: 93.0% CV -> omega^2 = log(1 + 0.93^2)

    # ============================================================
    # Inter-occasion variability on Vc/F and kel. Three occasions per
    # subject (the three-way crossover); the paper reports a single
    # IOV variance per parameter, so occasions 2 and 3 are fixed
    # equal to occasion 1 (NONMEM $OMEGA BLOCK(1) SAME pattern).
    # ============================================================
    etaiov_lvc_1 ~ 0.23036                                        # Table 3 IOV Vc/F: 50.9% CV -> omega^2 = log(1 + 0.509^2)
    etaiov_lvc_2 ~ fixed(0.23036)                                 # SAME as occasion 1
    etaiov_lvc_3 ~ fixed(0.23036)                                 # SAME as occasion 1
    etaiov_lkel_1 ~ 0.18365                                       # Table 3 IOV kel:  44.9% CV -> omega^2 = log(1 + 0.449^2)
    etaiov_lkel_2 ~ fixed(0.18365)                                # SAME as occasion 1
    etaiov_lkel_3 ~ fixed(0.18365)                                # SAME as occasion 1

    # ============================================================
    # Residual error. Back 2018 Table 3 reports a proportional error
    # of 0.608 (i.e. 60.8% CV) as the only residual error term.
    # ============================================================
    propSd <- 0.608; label("Proportional residual error (fraction)")     # Table 3 final: proportional = 0.608 (RSE 6.4%)
  })

  model({
    # ------------------------------------------------------------
    # Occasion indicators for IOV on Vc/F and kel. OCC = 1 fasted,
    # OCC = 2 standard meal, OCC = 3 high-fat meal in the Back 2018
    # three-way crossover.
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_lvc <- oc1 * etaiov_lvc_1 + oc2 * etaiov_lvc_2 + oc3 * etaiov_lvc_3
    iov_lkel <- oc1 * etaiov_lkel_1 + oc2 * etaiov_lkel_2 + oc3 * etaiov_lkel_3

    # ------------------------------------------------------------
    # Meal-type indicators. FED * (1 - FED_HIGHFAT) is 1 only for
    # the standard meal arm; FED_HIGHFAT is 1 only for the high-fat
    # meal arm; fasted records have both zero.
    # ------------------------------------------------------------
    fed_std <- FED * (1 - FED_HIGHFAT)
    fed_hf  <- FED_HIGHFAT

    # ------------------------------------------------------------
    # Time-varying gastric-emptying boost. Active during the first
    # mtime2 = 6.94 h after the drug dose (drug and meal coincide
    # in the Back 2018 study design); gated by FED so fasted
    # occasions never receive the boost.
    # ------------------------------------------------------------
    food_kg_window <- FED * ((t - tlast(stomach)) < mtime2)
    kg <- exp(lkg + etalkg) * (1 + e_food_kg * food_kg_window)

    # ------------------------------------------------------------
    # Other individual parameters. kg' shares the kg IIV per the
    # paper (Table 3 footnote "k_g (k_g', shared IIV)"); km&a, kout
    # are typical-value-only.
    # ------------------------------------------------------------
    kgp  <- exp(lkgp + etalkg)
    kma  <- exp(lkma)
    kel  <- exp(lkel + etalkel + iov_lkel)
    kout <- exp(lkout)
    vc   <- exp(lvc + etalvc + iov_lvc) * (1 + e_fed_vc_std * fed_std + e_fed_vc_hf * fed_hf)

    # ------------------------------------------------------------
    # Drug GI / disposition ODE system (Back 2018 Eqs. 1-3).
    # Stomach -> duodenum -> central; central elimination is
    # first-order; the duodenum -> central transfer is enhanced by
    # bile-acid (calorie-proxied) co-secretion through the
    # multiplicative term (1 + e_cal_kma * duodenum_food).
    # ------------------------------------------------------------
    kma_eff <- kma * (1 + e_cal_kma * duodenum_food)
    d/dt(stomach)  <- -kg * stomach
    d/dt(duodenum) <-  kg * stomach - kma_eff * duodenum
    d/dt(central)  <-  kma_eff * duodenum - kel * central

    # ------------------------------------------------------------
    # Calorie sub-model (Back 2018 Methods 'Pharmacokinetic
    # modeling', calorie ODEs). Stomach calories empty into the
    # duodenum with kgp = kg'; duodenum calories are eliminated
    # first-order with kout. Calories enter the system through a
    # dose into stomach_food whose amount is the kcal load.
    # ------------------------------------------------------------
    d/dt(stomach_food)  <- -kgp * stomach_food
    d/dt(duodenum_food) <-  kgp * stomach_food - kout * duodenum_food

    # ------------------------------------------------------------
    # Observation. central is in mg, vc is in L -> central / vc is
    # mg/L which equals ug/mL (Back 2018 reports Cmax and AUC in
    # ug/mL and ug.h/mL respectively; Table 4).
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
