Cheng_2026_levamisole_2cm <- function() {
  description <- paste(
    "Preclinical + human (duck, rabbit, chicken, goat, dog, sheep, pig, human).",
    "Interspecies meta-analysis two-compartment PK model for levamisole, jointly fitted to",
    "digitized intravenous plasma profiles from 8 species plus human oral data. The",
    "distribution parameters are shared across species through simple allometric relationships",
    "on body weight (Vp = a1*BW^b1, Vt = a2*BW^b2, CLd = a3*BW^b3), while elimination clearance",
    "is species-specific and selected by mutually exclusive SPECIES_* indicators, with human as",
    "the reference species. Absorption is first order into a depot; only the human oral arm was",
    "informed by the joint fit, so the depot ka and bioavailability are the human oral values.",
    "Naive-pooled deterministic fit in ADAPT with no between-subject variability.",
    "Companion model to Cheng_2026_levamisole_mpbpk.",
    sep = " "
  )
  reference <- paste(
    "Cheng C, Jeong YS, Jusko WJ.",
    "Meta-analysis of levamisole absorption and disposition across diverse species using a",
    "minimal physiologically-based pharmacokinetic model.",
    "J Pharm Investig. 2026;56(2):171-183.",
    "doi:10.1007/s40005-025-00770-6.",
    sep = " "
  )
  vignette <- "Cheng_2026_levamisole_interspecies"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    weight = "kg"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Levamisole plasma concentrations were the only measured
  # quantity in every source study pooled by Cheng 2026, so the disposition
  # states are plasma and the depot is the oral administration site.
  compartmentData <- list(
    depot       = list(analyte = "levamisole", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levamisole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levamisole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives every distribution parameter through the paper's allometric relations",
        "(Fig. 2 and Table 2): Vp = 0.339*WT^1.309 L, Vt = 0.937*WT^0.885 L and",
        "CLd = 0.352*WT^0.920 L/h. Because the intercepts a1, a2 and a3 are the values at",
        "WT = 1 kg, the model uses WT directly rather than a normalised WT/wtRef ratio.",
        "Also converts the weight-normalised species clearances of Table 1 (L/h/kg) into L/h.",
        "Body weights used by the joint fit were one typical value per species",
        "(Table 1 and Supplemental Table S1C): duck 2.5, rabbit 3, chicken 4.5, goat 18,",
        "dog 20.7, sheep 26, pig 39.2 and human 70 kg.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    SPECIES_DUCK = list(
      description        = "Duck species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = duck, 0 = another species in the model. Selects the duck clearance of",
        "Table 1 (joint 2CM column). Mutually exclusive with the other SPECIES_* indicators;",
        "all of them 0 selects the human parameter set. Ducks were the species Cheng 2026",
        "excluded from the simple allometric clearance regression as an outlier, which is why",
        "clearance is carried per species rather than allometrically scaled in this model.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_RABBIT = list(
      description        = "Rabbit species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = rabbit, 0 = another species in the model. Selects the rabbit clearance of",
        "Table 1 (joint 2CM column). Cheng 2026 does not report the rabbit strain.",
        "Mutually exclusive with the other SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_CHICKEN = list(
      description        = "Chicken species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = chicken, 0 = another species in the model. Selects the chicken clearance of",
        "Table 1 (joint 2CM column). The source chicken data are broiler breeder hens",
        "(Supplemental Table S2, El-Kholy 2006). Mutually exclusive with the other",
        "SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_GOAT = list(
      description        = "Goat species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = goat, 0 = another species in the model. Selects the goat clearance of",
        "Table 1 (joint 2CM column). Mutually exclusive with the other SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_DOG = list(
      description        = "Dog species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = dog, 0 = another species in the model. Selects the dog clearance of",
        "Table 1 (joint 2CM column). The dog intravenous profile declines essentially",
        "mono-exponentially (Supplemental Figure S1), so the shared allometric distribution",
        "parameters are poorly informed by this species. Mutually exclusive with the other",
        "SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_SHEEP = list(
      description        = "Sheep species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = sheep, 0 = another species in the model. Selects the sheep clearance of",
        "Table 1 (joint 2CM column). Mutually exclusive with the other SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_PIG = list(
      description        = "Pig species indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pig, 0 = another species in the model. Selects the pig clearance of",
        "Table 1 (joint 2CM column). Cheng 2026 reports that the pig fit shows notable",
        "discrepancies under the joint allometric 2CM (Results, Two-compartment model",
        "fitting). Mutually exclusive with the other SPECIES_* indicators.",
        sep = " "
      ),
      source_name        = "SPECIES"
    )
  )

  population <- list(
    species        = paste(
      "duck + rabbit + chicken + goat + dog + sheep + pig + human",
      "(joint interspecies fit; human is the model's reference species)",
      sep = " "
    ),
    n_subjects     = NA_integer_,
    n_studies      = 8L,
    weight_range   = "2.5-70 kg (one typical body weight per species; Table 1)",
    disease_state  = "Healthy / experimentally infected animals and healthy human volunteers",
    dose_range     = paste(
      "Intravenous single doses of 5-40 mg/kg in the animal species and single oral doses of",
      "50-150 mg in humans (Supplemental Tables S2, S3 and S4).",
      sep = " "
    ),
    notes          = paste(
      "Meta-analysis of published levamisole PK. Concentration-time profiles were digitized",
      "from over 40 publications covering 18 species with WebPlotDigitizer; 8 species met the",
      "inclusion criteria for joint modelling (availability of intravenous data except for",
      "humans, and exclusion of cold-blooded species). One representative single-dose profile",
      "was used per species so that each contributes equal weight, with the middle dose chosen",
      "for species studied over a dose range (rabbits, sheep). The fit is naive-pooled and",
      "deterministic (maximum likelihood in ADAPT 5 on the pooled typical profiles), so there",
      "is no between-subject variability to report; n_subjects is therefore not applicable.",
      "The source studies for each species are listed in Supplemental Tables S2 (intravenous),",
      "S3 (extravascular) and S4 (human oral).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Allometric distribution parameters (Fig. 2, Table 2). Cheng 2026 writes
    # these as Y = a * BW^b with body weight in kg, so each intercept a IS the
    # parameter value at 1 kg and is stored here as the log-scale typical
    # value with the exponent carried alongside.
    #
    # Naming map from the paper to the library canonicals: the paper's central
    # volume Vp is `vc`, its peripheral volume Vt is `vp`, and its
    # distributional clearance CLd is `q`.
    # ---------------------------------------------------------------------
    lvc     <- log(0.339) ; label("Central volume Vp at 1 kg body weight (L)")                # Table 2 a1 = 0.339 (CV% 12.68)
    e_wt_vc <- 1.309      ; label("Allometric exponent on the central volume Vp (unitless)")  # Table 2 b1 = 1.309 (CV% 3.33)
    lvp     <- log(0.937) ; label("Peripheral volume Vt at 1 kg body weight (L)")             # Table 2 a2 = 0.937 (CV% 21.13)
    e_wt_vp <- 0.885      ; label("Allometric exponent on the peripheral volume Vt (unitless)") # Table 2 b2 = 0.885 (CV% 8.33)
    lq      <- log(0.352) ; label("Distributional clearance CLd at 1 kg body weight (L/h)")   # Table 2 a3 = 0.352 (CV% 28.20)
    e_wt_q  <- 0.920      ; label("Allometric exponent on the distributional clearance CLd (unitless)") # Table 2 b3 = 0.920 (CV% 10.89)

    # ---------------------------------------------------------------------
    # Species-specific elimination clearance, Table 1 "Joint 2CM" column, in
    # L/h/kg; model() multiplies by WT to obtain L/h. Cheng 2026 deliberately
    # did NOT scale clearance allometrically in the joint fits: the reported
    # exponent b = 0.26 had too low an R-squared (0.46, Fig. 4) to be usable.
    # Human is the reference species (all SPECIES_* indicators 0).
    # ---------------------------------------------------------------------
    lcl         <- log(0.25)  ; label("Elimination clearance, human (L/h/kg)")   # Table 1 human, joint 2CM CL = 0.25 (CV% 3.72)
    lcl_duck    <- log(0.20)  ; label("Elimination clearance, duck (L/h/kg)")    # Table 1 duck, joint 2CM CL = 0.20 (CV% 7.54)
    lcl_rabbit  <- log(1.09)  ; label("Elimination clearance, rabbit (L/h/kg)")  # Table 1 rabbit, joint 2CM CL = 1.09 (CV% 7.79)
    lcl_chicken <- log(0.93)  ; label("Elimination clearance, chicken (L/h/kg)") # Table 1 chicken, joint 2CM CL = 0.93 (CV% 9.07)
    lcl_goat    <- log(0.38)  ; label("Elimination clearance, goat (L/h/kg)")    # Table 1 goat, joint 2CM CL = 0.38 (CV% 0.90)
    lcl_dog     <- log(0.49)  ; label("Elimination clearance, dog (L/h/kg)")     # Table 1 dog, joint 2CM CL = 0.49 (CV% 4.39)
    lcl_sheep   <- log(0.99)  ; label("Elimination clearance, sheep (L/h/kg)")   # Table 1 sheep, joint 2CM CL = 0.99 (CV% 5.88)
    lcl_pig     <- log(0.27)  ; label("Elimination clearance, pig (L/h/kg)")     # Table 1 pig, joint 2CM CL = 0.27 (CV% 7.43)

    # ---------------------------------------------------------------------
    # Absorption (Eq. 12). Humans contributed no intravenous data, so the
    # joint fit estimated the human oral ka and fixed the human oral F at the
    # literature value; these are therefore HUMAN oral values, not shared
    # cross-species ones. Per-species oral ka and F from the separate
    # individual-species fits are in Supplemental Table S6 and are not part of
    # this joint model; see the vignette Errata.
    # ---------------------------------------------------------------------
    lka     <- log(1.702)        ; label("First-order oral absorption rate, human (1/h)")  # Table 2 ka = 1.702 (CV% 12.10)
    lfdepot <- fixed(log(0.66))  ; label("Oral bioavailability, human (fraction)")         # Table 2 F = 0.66 (Fixed), literature value; Supplemental Table S9 reports 62.5-68% in humans

    # ---------------------------------------------------------------------
    # Residual error. Cheng 2026 fitted with the ADAPT variance model
    # Var(i) = (sigma1 + sigma2 * Y(i))^2, i.e. a combined additive plus
    # proportional error on the standard-deviation scale (Methods, Model
    # fitting), but reports neither sigma1 nor sigma2 anywhere in the paper or
    # its supplement. The FORM is therefore preserved and both magnitudes are
    # fixed at zero rather than invented; see the vignette Errata.
    # ---------------------------------------------------------------------
    addSd  <- fixed(0) ; label("Additive residual error, sigma1 (mg/L; not published by Cheng 2026)")
    propSd <- fixed(0) ; label("Proportional residual error, sigma2 (fraction; not published by Cheng 2026)")
  })

  model({
    # -------------------------------------------------------------------
    # 1. Species selection. The indicators are mutually exclusive, so the
    #    human (reference) flag is one minus their sum. A data set that sets
    #    two indicators to 1 makes isHuman negative and produces nonsense
    #    rather than erroring -- assemblers must preserve exclusivity.
    # -------------------------------------------------------------------
    isDuck    <- SPECIES_DUCK
    isRabbit  <- SPECIES_RABBIT
    isChicken <- SPECIES_CHICKEN
    isGoat    <- SPECIES_GOAT
    isDog     <- SPECIES_DOG
    isSheep   <- SPECIES_SHEEP
    isPig     <- SPECIES_PIG
    isHuman   <- 1 - isDuck - isRabbit - isChicken - isGoat - isDog - isSheep - isPig

    # -------------------------------------------------------------------
    # 2. Individual parameters. Distribution is allometric and shared;
    #    elimination clearance is species-specific and weight-normalised.
    # -------------------------------------------------------------------
    vc <- exp(lvc) * WT^e_wt_vc   # Table 2: Vp = a1 * BW^b1
    vp <- exp(lvp) * WT^e_wt_vp   # Table 2: Vt = a2 * BW^b2
    q  <- exp(lq)  * WT^e_wt_q    # Table 2: CLd = a3 * BW^b3

    clPerKg <-
      isDuck    * exp(lcl_duck) +
      isRabbit  * exp(lcl_rabbit) +
      isChicken * exp(lcl_chicken) +
      isGoat    * exp(lcl_goat) +
      isDog     * exp(lcl_dog) +
      isSheep   * exp(lcl_sheep) +
      isPig     * exp(lcl_pig) +
      isHuman   * exp(lcl)
    cl <- clPerKg * WT            # Table 1 reports CL in L/h/kg

    ka <- exp(lka)

    # -------------------------------------------------------------------
    # 3. ODE system, Eqs. (3), (4) and (12). Cheng 2026 writes the
    #    disposition equations on the concentration scale
    #    (Vp*dCp/dt = CLd*(Ct - Cp) - CL*Cp; Vt*dCt/dt = CLd*(Cp - Ct));
    #    the amount-scale form below is algebraically identical, with the
    #    absorption rate F*ka*Aa entering the central compartment via the
    #    depot and f(depot).
    # -------------------------------------------------------------------
    Cc <- central / vc            # paper's Cp
    Ct <- peripheral1 / vp        # paper's Ct

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + q * (Ct - Cc) - cl * Cc
    d/dt(peripheral1) <-  q * (Cc - Ct)

    f(depot) <- exp(lfdepot)

    Cc ~ add(addSd) + prop(propSd)
  })
}
