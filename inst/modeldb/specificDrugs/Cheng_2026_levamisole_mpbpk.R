Cheng_2026_levamisole_mpbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based, interspecies).",
    "Preclinical + human (duck, rabbit, chicken, goat, dog, sheep, pig, human).",
    "Levamisole minimal PBPK model jointly fitted to digitized intravenous plasma profiles from",
    "8 species plus human oral data. A blood compartment exchanges with two lumped tissue pools",
    "that split cardiac output by the fractional distribution parameters fd1 and fd2 and split",
    "the residual body volume by the tissue-volume fraction ft. Distribution is perfusion",
    "limited (fd1 + fd2 = 1) and governed by a single tissue-to-plasma partition coefficient Kp",
    "shared across species, with species-specific Kp values for pig and chicken. Cardiac output",
    "is allometric (Qco = 14.1*BW^0.75 L/h, Brown 1997) and blood volume is a species-specific",
    "fraction of body weight; elimination clearance is species-specific and selected by mutually",
    "exclusive SPECIES_* indicators, with human as the reference species. Absorption is first",
    "order into a depot and was informed only by the human oral arm. Naive-pooled deterministic",
    "fit in ADAPT with no between-subject variability. Cheng 2026 preferred this model over the",
    "allometric two-compartment companion Cheng_2026_levamisole_2cm on AIC.",
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
  # biological matrix. `central` is Cheng 2026's Blood compartment (Fig. 3):
  # its amount divided by Vb*Rb is the PLASMA concentration Cp, which is the
  # quantity every pooled source study measured. `peripheral1` and
  # `peripheral2` are the paper's lumped Tissue 1 and Tissue 2 pools; they are
  # anatomical tissue volumes rather than a single named organ, because Cheng
  # 2026 explicitly could not apply a tissue-lumping assignment for levamisole
  # (Discussion, Distribution kinetics: no tissue-specific fd or Kp data exist).
  compartmentData <- list(
    depot       = list(analyte = "levamisole", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levamisole", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "levamisole", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "levamisole", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the whole physiology of the model. Cardiac output Qco = 14.1*WT^0.75 L/h",
        "(Eq. 8, Brown 1997; reproduced per species in Supplemental Table S1C). Blood volume",
        "Vb = fvb*WT with a species-specific fraction fvb from Supplemental Table S1A. Assuming",
        "unit tissue density, the two lumped tissue volumes are V1 = (WT - Vb)*ft and",
        "V2 = (WT - Vb)*(1 - ft) (Eqs. 9 and 10), so WT in kg maps directly onto volumes in L.",
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
        "1 = duck, 0 = another species in the model. Selects the duck clearance of Table 1",
        "(joint mPBPK column) and the duck blood-volume fraction of 86.3 mL/kg reported",
        "directly in Supplemental Table S1A. Ducks take the shared Kp of 1.41, not a",
        "species-specific one. Mutually exclusive with the other SPECIES_* indicators;",
        "all of them 0 selects the human parameter set.",
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
        "Table 1 (joint mPBPK column) and the rabbit blood-volume fraction, which is derived",
        "as the 152.7 mL reported in Supplemental Table S1A divided by the paper's 3 kg rabbit",
        "body weight. Rabbits take the shared Kp of 1.41. Cheng 2026 does not report the",
        "rabbit strain. Mutually exclusive with the other SPECIES_* indicators.",
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
        "Table 1 (joint mPBPK column), the chicken blood-volume fraction of 10 percent of body",
        "weight from Supplemental Table S1A, AND the chicken-specific tissue-to-plasma",
        "partition coefficient Kp,chicken = 3.38 from Table 3. Cheng 2026 needed the separate",
        "Kp because the chicken steady-state volume of distribution (8.07 L/kg by NCA) is far",
        "above the cross-species norm despite unremarkable plasma protein binding (19.4%).",
        "Mutually exclusive with the other SPECIES_* indicators.",
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
        "1 = goat, 0 = another species in the model. Selects the goat clearance of Table 1",
        "(joint mPBPK column) and the goat blood-volume fraction of 70 mL/kg from Supplemental",
        "Table S1A. Goats take the shared Kp of 1.41. Mutually exclusive with the other",
        "SPECIES_* indicators.",
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
        "1 = dog, 0 = another species in the model. Selects the dog clearance of Table 1",
        "(joint mPBPK column) and the dog blood-volume fraction of 84 mL/kg from Supplemental",
        "Table S1A. Dogs take the shared Kp of 1.41. The dog intravenous profile declines",
        "essentially mono-exponentially, which Cheng 2026 identifies as the reason for the poor",
        "CV% on the dog distribution parameters (Discussion, Model comparisons and performance).",
        "Mutually exclusive with the other SPECIES_* indicators.",
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
        "1 = sheep, 0 = another species in the model. Selects the sheep clearance of Table 1",
        "(joint mPBPK column) and the sheep blood-volume fraction of 59 mL/kg from Supplemental",
        "Table S1A. Sheep take the shared Kp of 1.41. Mutually exclusive with the other",
        "SPECIES_* indicators.",
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
        "1 = pig, 0 = another species in the model. Selects the pig clearance of Table 1",
        "(joint mPBPK column), the pig blood-volume fraction of 60 mL/kg from Supplemental",
        "Table S1A, AND the pig-specific tissue-to-plasma partition coefficient Kp,pig = 5.62",
        "from Table 3. As for chickens, the separate Kp reflects a steady-state volume of",
        "distribution (4.73 L/kg by NCA) that plasma protein binding alone does not explain.",
        "Mutually exclusive with the other SPECIES_* indicators.",
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
      "humans, and exclusion of cold-blooded species, whose basal metabolic rate is much",
      "lower). One representative single-dose profile was used per species so that each",
      "contributes equal weight. The fit is naive-pooled and deterministic (maximum likelihood",
      "in ADAPT 5 on the pooled typical profiles), so there is no between-subject variability",
      "to report; n_subjects is therefore not applicable. This model gave a lower AIC",
      "(598,377) than the allometric two-compartment companion model (598,512).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Distribution parameters (Table 3). Perfusion-limited: the assumption of
    # low tissue permeability gave unsuitable fits, so fd1 + fd2 was fixed at
    # its upper bound of 1 (Eq. 11, Results, Minimal PBPK model fitting).
    # ---------------------------------------------------------------------
    lkp             <- log(1.41) ; label("Tissue-to-plasma partition coefficient, all species except pig and chicken (unitless)") # Table 3 Kp = 1.41 (CV% 2.40)
    lkp_pig         <- log(5.62) ; label("Tissue-to-plasma partition coefficient, pig (unitless)")     # Table 3 Kp,pig = 5.62 (CV% 7.71); the Results text prints the rounded 5.61, see vignette Errata
    lkp_chicken     <- log(3.38) ; label("Tissue-to-plasma partition coefficient, chicken (unitless)") # Table 3 Kp,chicken = 3.38 (CV% 4.47)
    frac_tissue1    <- 0.497     ; label("Fraction of the total tissue volume assigned to tissue 1, ft (unitless)") # Table 3 ft = 0.497 (CV% 2.71)
    frac_qco_tissue1 <- 0.038    ; label("Fraction of cardiac output distributed to tissue 1, fd1 (unitless)")      # Table 3 fd1 = 0.038 (CV% 5.75)
    fd_total        <- fixed(1)  ; label("Sum of the fractional distribution parameters, fd,total (unitless)")      # Table 3 fd,total = 1.0; Eq. 11 upper bound, perfusion-limited distribution. fd2 = fd_total - fd1 = 0.962 (Table 3)

    # ---------------------------------------------------------------------
    # Physiological constants.
    # ---------------------------------------------------------------------
    lrb      <- fixed(log(1))  ; label("Blood-to-plasma partition coefficient, Rb (unitless)")    # Methods, Calculation of basic PK parameters: "Rb ... was assumed to be 1"
    qco_ref  <- fixed(14.1)    ; label("Cardiac output at 1 kg body weight (L/h)")                # Eq. 8: Qco (L/h) = 14.1 * BW(kg)^0.75 (Brown 1997)
    e_wt_qco <- fixed(0.75)    ; label("Allometric exponent on cardiac output (unitless)")        # Eq. 8 exponent 0.75

    # ---------------------------------------------------------------------
    # Species blood volume as a fraction of body weight, from the literature
    # values compiled in Supplemental Table S1A. Six species are reported
    # directly as a fraction (mL/kg, or 10% of body weight for the chicken);
    # rabbit and human are reported as absolute volumes and are converted here
    # by dividing by the body weight the joint fit used for that species, so
    # each value reproduces the paper's Vb exactly at the paper's body weight.
    # ---------------------------------------------------------------------
    fvb         <- fixed(0.0714) ; label("Blood volume fraction, human (L/kg)")   # Supplemental Table S1A human 5 L / 70 kg = 0.0714 L/kg (Sharma & Sharma 2018)
    fvb_duck    <- fixed(0.0863) ; label("Blood volume fraction, duck (L/kg)")    # Supplemental Table S1A duck 86.3 mL/kg (Portman 1952)
    fvb_rabbit  <- fixed(0.0509) ; label("Blood volume fraction, rabbit (L/kg)")  # Supplemental Table S1A rabbit 152.7 mL / 3 kg = 0.0509 L/kg (Baby 2014)
    fvb_chicken <- fixed(0.1)    ; label("Blood volume fraction, chicken (L/kg)") # Supplemental Table S1A chicken "10%" of body weight (Newell & Shaffner 1950)
    fvb_goat    <- fixed(0.070)  ; label("Blood volume fraction, goat (L/kg)")    # Supplemental Table S1A goat 70 mL/kg (Courtice 1943)
    fvb_dog     <- fixed(0.084)  ; label("Blood volume fraction, dog (L/kg)")     # Supplemental Table S1A dog 84 mL/kg (Gibson 1938)
    fvb_sheep   <- fixed(0.059)  ; label("Blood volume fraction, sheep (L/kg)")   # Supplemental Table S1A sheep 59 mL/kg (Hansard 1956)
    fvb_pig     <- fixed(0.060)  ; label("Blood volume fraction, pig (L/kg)")     # Supplemental Table S1A pig 60 mL/kg (Bush 1955)

    # ---------------------------------------------------------------------
    # Species-specific elimination clearance, Table 1 "Joint mPBPK" column, in
    # L/h/kg; model() multiplies by WT to obtain L/h. Human is the reference
    # species (all SPECIES_* indicators 0).
    # ---------------------------------------------------------------------
    lcl         <- log(0.249)  ; label("Elimination clearance, human (L/h/kg)")   # Table 1 human, joint mPBPK CL = 0.249 (CV% 5.52)
    lcl_duck    <- log(0.216)  ; label("Elimination clearance, duck (L/h/kg)")    # Table 1 duck, joint mPBPK CL = 0.216 (CV% 4.25)
    lcl_rabbit  <- log(1.46)   ; label("Elimination clearance, rabbit (L/h/kg)")  # Table 1 rabbit, joint mPBPK CL = 1.46 (CV% 6.37)
    lcl_chicken <- log(2.59)   ; label("Elimination clearance, chicken (L/h/kg)") # Table 1 chicken, joint mPBPK CL = 2.59 (CV% 4.05)
    lcl_goat    <- log(0.359)  ; label("Elimination clearance, goat (L/h/kg)")    # Table 1 goat, joint mPBPK CL = 0.359 (CV% 1.05)
    lcl_dog     <- log(0.464)  ; label("Elimination clearance, dog (L/h/kg)")     # Table 1 dog, joint mPBPK CL = 0.464 (CV% 4.12)
    lcl_sheep   <- log(0.994)  ; label("Elimination clearance, sheep (L/h/kg)")   # Table 1 sheep, joint mPBPK CL = 0.994 (CV% 4.56)
    lcl_pig     <- log(0.371)  ; label("Elimination clearance, pig (L/h/kg)")     # Table 1 pig, joint mPBPK CL = 0.371 (CV% 7.28)

    # ---------------------------------------------------------------------
    # Absorption (Eq. 12), human oral arm only. Table 3 does NOT report the
    # absorption rate constant of the joint mPBPK fit, so the value used here
    # is the human ka of the individual-species mPBPK fit reported in the
    # footnote to Supplemental Table S7. Bioavailability was fixed at 0.66 in
    # both fits. See the vignette Errata.
    # ---------------------------------------------------------------------
    lka     <- log(1.156)       ; label("First-order oral absorption rate, human (1/h)")  # Supplemental Table S7 footnote b: human ka fitted to 1.156 /h with F fixed at 0.66
    lfdepot <- fixed(log(0.66)) ; label("Oral bioavailability, human (fraction)")         # Absorption kinetics: human F fixed at 0.66 from the literature; Supplemental Table S9 reports 62.5-68% in humans

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
    # 2. Species-selected parameters. Only pig and chicken carry their own
    #    partition coefficient; every other species uses the shared Kp.
    # -------------------------------------------------------------------
    kp <-
      isPig     * exp(lkp_pig) +
      isChicken * exp(lkp_chicken) +
      (1 - isPig - isChicken) * exp(lkp)

    fvbSpecies <-
      isDuck    * fvb_duck +
      isRabbit  * fvb_rabbit +
      isChicken * fvb_chicken +
      isGoat    * fvb_goat +
      isDog     * fvb_dog +
      isSheep   * fvb_sheep +
      isPig     * fvb_pig +
      isHuman   * fvb

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

    # -------------------------------------------------------------------
    # 3. Physiology. Eq. (8) for cardiac output, Eqs. (9) and (10) for the
    #    two lumped tissue volumes assuming unit tissue density, and Eq. (11)
    #    for the split of cardiac output between them.
    # -------------------------------------------------------------------
    rb   <- exp(lrb)
    qco  <- qco_ref * WT^e_wt_qco
    vb   <- fvbSpecies * WT
    ft   <- frac_tissue1
    fd1  <- frac_qco_tissue1
    fd2  <- fd_total - fd1
    v1   <- (WT - vb) * ft
    v2   <- (WT - vb) * (1 - ft)

    ka <- exp(lka)

    # -------------------------------------------------------------------
    # 4. ODE system, Eqs. (5), (6), (7) and (12). Cheng 2026 writes these on
    #    the concentration scale with the blood-compartment left-hand side
    #    Vb*Rb*dCp/dt; the amount-scale form below is algebraically identical.
    #    Dividing the `central` amount by Vb*Rb therefore recovers the PLASMA
    #    concentration Cp, and an intravenous bolus into `central` reproduces
    #    the paper's initial condition Cp(0) = Dose/(Vb*Rb) exactly.
    # -------------------------------------------------------------------
    Cc <- central / (vb * rb)     # paper's Cp, the plasma concentration
    C1 <- peripheral1 / v1        # paper's C1, tissue 1 concentration
    C2 <- peripheral2 / v2        # paper's C2, tissue 2 concentration

    d/dt(depot) <- -ka * depot

    d/dt(central) <- ka * depot +
      qco * fd1 * rb * (C1 / kp - Cc) +
      qco * fd2 * rb * (C2 / kp - Cc) -
      cl * Cc

    d/dt(peripheral1) <- qco * fd1 * rb * (Cc - C1 / kp)
    d/dt(peripheral2) <- qco * fd2 * rb * (Cc - C2 / kp)

    f(depot) <- exp(lfdepot)

    Cc ~ add(addSd) + prop(propSd)
  })
}
