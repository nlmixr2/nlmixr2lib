Mulford_2026_vonoprazan_maturation <- function() {
  description <- paste0(
    "Organ-maturation extension of the Mulford 2026 Model 3 vonoprazan ",
    "population PK model, used to project steady-state exposure in infants ",
    "and young children aged 1 month to under 6 years. The structural, ",
    "covariate, variability and residual-error model is identical to ",
    "modellib('Mulford_2026_vonoprazan') -- two compartments, three ",
    "transit compartments at 4/MTT followed by first-order absorption, a ",
    "dose-power and female effect on relative bioavailability, female, ",
    "age-power and day-after-first-dose effects on absorption rate, a ",
    "body-weight power effect on central volume and day-after-first-dose ",
    "effects on clearance and central volume -- with one addition: ",
    "clearance is multiplied by a sigmoid Hill maturation factor of ",
    "postmenstrual age, CLmat = Fmat * CL with ",
    "Fmat = pma^3.4 / (47.3^3.4 + pma^3.4), where postmenstrual age in ",
    "weeks is derived from chronological age in years as ",
    "pma = age * 365.25/7 + 42. The 47.3-week half-maturation age and the ",
    "Hill coefficient of 3.4 are literature constants, not estimated by ",
    "Mulford 2026; no observed vonoprazan data below 6 years of age exist. ",
    "Fmat is essentially 1 by school age, so this model reduces to the ",
    "parent model in adolescents and adults. Simulations from it produced ",
    "Mulford 2026 Figure 6, which supports 10 mg and 20 mg doses in ",
    "children aged 1 to under 6 years and suggests lower doses for infants ",
    "aged 1 month to under 1 year."
  )
  reference <- paste(
    "Mulford DJ, Facius A, Witt G, Howden CW, Wagner T, Leifke E,",
    "Scarpignato C. Development and Use of a Population Pharmacokinetic",
    "Model for Characterizing the Pharmacokinetics of Vonoprazan in",
    "Pediatric Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70291.",
    "doi:10.1002/psp4.70291.",
    "Maturation extension from Methods 2.4.2; structural detail and",
    "covariate centering values from the Supporting Information Model Code",
    "(the final-model NONMEM control stream, PSP4-15-e70291-s001.docx).",
    "Parent model: modellib('Mulford_2026_vonoprazan').",
    sep = " "
  )
  vignette <- "Mulford_2026_vonoprazan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model on central volume only, centered at 70 kg (Mulford",
        "2026 Methods 2.2.2; Supporting Information Model Code line",
        "TVVC = TVVC * (WEIGHT/70)**THETA(11)). For the Figure 6",
        "simulations, Methods 2.4.2 sampled age- and sex-matched body",
        "weights from the WHO/CDC normal growth statistics. No weight",
        "effect on clearance was retained (Results 3.3.3)."
      ),
      source_name = "WEIGHT"
    ),
    AGE = list(
      description = "Chronological age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used TWICE in this model, and the two uses have different",
        "centerings and different mechanisms. (1) As a power model on the",
        "absorption rate constant centered at 28 years, per the Supporting",
        "Information Model Code line TVKA = TVKA * (AGE/28)**THETA(14).",
        "Mulford 2026 Methods 2.2.2 states 18 years for this centering,",
        "which contradicts the deposited control stream; the control-stream",
        "value is used here. See the vignette Errata. (2) As the input to",
        "the organ-maturation factor on clearance, converted to",
        "postmenstrual age in weeks inside model() as",
        "pma = AGE * 365.25/7 + 42 (Methods 2.4.2). The conversion assumes",
        "a 42-week postmenstrual age at birth, so a term-born infant of",
        "chronological age 0 starts at pma = 42 weeks; a preterm infant is",
        "not distinguished. This model is calibrated on observed data only",
        "down to 6 years of age -- everything below that is extrapolation",
        "carried entirely by the maturation function."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The source column is SEX with SEX == 2 denoting female",
        "(Supporting Information Model Code: IF( SEX == 2 )), so",
        "SEXF = as.integer(SEX == 2). The Methods 2.4.2 virtual population",
        "was 50 percent boys and 50 percent girls."
      ),
      source_name = "SEX"
    ),
    DAY2 = list(
      description = "Day-after-first-dose landmark indicator: 1 = the record falls on study day 2 or later, 0 = the record falls on study day 1 (the first dosing day).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study day 1, the first dosing day)",
      notes = paste(
        "Implements the Supporting Information Model Code branches",
        "IF(DAY>1) on absorption rate, clearance and central volume, i.e.",
        "DAY2 = as.integer(study_day >= 2). Because the Figure 6",
        "projections are steady-state exposures, DAY2 = 1 is the relevant",
        "setting for reproducing them."
      ),
      source_name = "DAY"
    ),
    DOSE_VONOPRAZAN_MG = list(
      description = "Administered vonoprazan dose level carried on every record of the dosing interval, in mg.",
      units = "mg",
      type = "continuous",
      reference_category = "20 mg (the dose at which relative bioavailability equals its typical value)",
      notes = paste(
        "Continuous power-model regressor on relative bioavailability,",
        "centered at 20 mg (Supporting Information Model Code line",
        "TVFREL = TVFREL * (DOSE/20)**THETA(10)). Methods 2.4.3 swept",
        "candidate pediatric doses from 1 mg to 25 mg through this term",
        "when matching pediatric to adult steady-state AUC; the source",
        "studies span 1 mg to 120 mg (Table S1)."
      ),
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    RACE = list(
      description = "Self-reported race category (Asian, Black, Other, White).",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Carried in the analysis dataset and screened graphically in Figure",
        "S3, but not retained; no point estimate is reported."
      ),
      source_name = "RACE"
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate.",
      units = "mL/min/1.73m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Present in the analysis dataset but not tested or retained in the",
        "published covariate model; no point estimate is reported."
      ),
      source_name = "EGFR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 392L,
    n_studies = 10L,
    n_observations = 8201L,
    age_range = "6-54 years observed; applied by extrapolation from 1 month of age",
    weight_range = "20.7-132 kg observed",
    weight_median = "66 kg (adults), 64 kg (adolescents), 32 kg (children)",
    disease_state = "healthy volunteers (all eight adult studies) and patients with gastroesophageal reflux disease (the adolescent study VPED-102 and the child study VPED-103)",
    dose_range = "1-120 mg oral observed; 1-25 mg swept in the pediatric dose-matching simulations",
    regions = "Japan, Europe, China, United States",
    notes = paste(
      "The estimation dataset is identical to that of",
      "modellib('Mulford_2026_vonoprazan') -- Mulford 2026 Tables S1 to S3",
      "-- because the maturation factor was ADDED to the fitted Model 3",
      "for simulation rather than estimated from data. The simulated",
      "population of Methods 2.4.2 was 50 percent boys and 50 percent",
      "girls aged 1 month to under 6 years in 1-month steps, 100000",
      "subjects per step for 14200000 virtual subjects in total, with age-",
      "and sex-matched body weights from the WHO/CDC normal growth",
      "statistics. No observed vonoprazan pharmacokinetic data below 6",
      "years of age exist, so every prediction from this model below that",
      "age rests on the assumed maturation function."
    )
  )

  ini({
    # Structural, covariate, variability and residual parameters are identical
    # to modellib('Mulford_2026_vonoprazan') and come from Mulford 2026
    # Table 1, rightmost (Children) column = Model 3. See that model file for
    # the per-parameter derivations; only the two maturation constants below
    # are new.
    lfdepot <- fixed(log(1)); label("Typical relative oral bioavailability at a 20 mg dose in males (unitless)")

    lmtt <- log(0.762); label("Typical mean transit time of the absorption delay chain (h)")           # Table 1 Children, MTT TV, RSE 3.3
    lka  <- log(3.08);  label("Typical first-order absorption rate constant from the third transit compartment (1/h)") # Table 1 Children, ka TV, RSE 7.7
    lcl  <- log(118);   label("Typical fully-matured apparent elimination clearance (L/h)")            # Table 1 Children, CL TV, RSE 2.2
    lvc  <- log(751);   label("Typical apparent central volume of distribution (L)")                   # Table 1 Children, Vc TV, RSE 2.8
    lq   <- log(49.8);  label("Typical apparent distribution clearance (L/h)")                         # Table 1 Children, Q TV, RSE 4.8
    lvp  <- log(271);   label("Typical apparent peripheral volume of distribution (L)")                # Table 1 Children, Vp TV, RSE 2.2

    e_dose_fdepot <- 0.290;  label("Power exponent on DOSE_VONOPRAZAN_MG/20 for relative bioavailability (unitless)") # Table 1 Children, Frel dose-effect, RSE 4.8
    e_sexf_fdepot <- 0.351;  label("Fractional change in relative bioavailability for females (unitless)")            # Table 1 Children, Frel female-effect 35.1, RSE 10.8
    e_sexf_ka     <- -0.490; label("Fractional change in absorption rate constant for females (unitless)")            # Table 1 Children, ka female-effect -49.0, RSE 12.9
    e_age_ka      <- -0.693; label("Power exponent on AGE/28 for the absorption rate constant (unitless)")            # Table 1 Children, ka age-effect, RSE 10.6
    e_day2_ka     <- -0.208; label("Fractional change in absorption rate constant on study day 2 and later (unitless)") # Table 1 Children, ka Day >1-effect -20.8, RSE 9.9
    e_day2_cl     <- -0.115; label("Fractional change in elimination clearance on study day 2 and later (unitless)")    # Table 1 Children, CL Day >1-effect -11.5, RSE 3.5
    e_wt_vc       <- 0.668;  label("Power exponent on WT/70 for the central volume of distribution (unitless)")       # Table 1 Children, Vc weight-effect, RSE 6.6
    e_day2_vc     <- -0.116; label("Fractional change in central volume on study day 2 and later (unitless)")         # Table 1 Children, Vc Day >1-effect -11.6, RSE 11.4

    # Organ-maturation constants. Mulford 2026 Methods 2.4.2 takes these from
    # a generic maturation function in the cited literature (references 22 and
    # 23) rather than estimating them, so both are held constant.
    pma_tm50 <- fixed(47.3); label("Postmenstrual age at which 50 percent of clearance maturation is achieved (weeks)") # Methods 2.4.2, printed maturation equation
    pma_hill <- fixed(3.4);  label("Hill coefficient describing the steepness of the clearance maturation function (unitless)") # Methods 2.4.2, printed maturation equation

    # Between-subject variability. See the parent model file for the argument
    # that the Table 1 BSV a rows are standard deviations. Variances are the
    # squares of the Table 1 Children BSV a entries; the covariance is
    # 0.913 * 0.366 * 0.389.
    etalcl + etalvc ~ c(0.133956,
                        0.129987, 0.151321)                                                            # Table 1 Children, CL BSV a 0.366, Cor 0.913, Vc BSV a 0.389
    etalmtt ~ 0.276676                                                                                 # Table 1 Children, MTT BSV a 0.526, RSE 4.1
    etalka  ~ 0.413449                                                                                 # Table 1 Children, ka BSV a 0.643, RSE 7.7

    propSd <- 0.247;        label("Proportional residual error (fraction)")   # Table 1 Children, Residual variability Prop 24.7, RSE 0.4
    addSd  <- fixed(0.001); label("Additive residual error (ng/mL)")          # Table 1 Children, Residual variability Add, reported as fixed
  })

  model({
    # 1. Derived covariate terms.
    dose_fdepot <- (DOSE_VONOPRAZAN_MG / 20)^e_dose_fdepot
    age_ka      <- (AGE / 28)^e_age_ka
    wt_vc       <- (WT / 70)^e_wt_vc

    #    Organ-maturation factor on clearance (Mulford 2026 Methods 2.4.2).
    #    Postmenstrual age in weeks is derived from chronological age in years
    #    assuming 42 weeks of gestation:  pma = age * 365.25/7 + 42.
    pma_wk        <- AGE * 365.25 / 7 + 42
    maturation_cl <- pma_wk^pma_hill / (pma_tm50^pma_hill + pma_wk^pma_hill)

    # 2. Individual parameters. CLmat = Fmat * CL (Methods 2.4.2).
    frel <- exp(lfdepot) * dose_fdepot * (1 + e_sexf_fdepot * SEXF)
    mtt  <- exp(lmtt + etalmtt)
    ka   <- exp(lka + etalka) * (1 + e_sexf_ka * SEXF) * age_ka * (1 + e_day2_ka * DAY2)
    cl   <- exp(lcl + etalcl) * (1 + e_day2_cl * DAY2) * maturation_cl
    vc   <- exp(lvc + etalvc) * wt_vc * (1 + e_day2_vc * DAY2)
    q    <- exp(lq)
    vp   <- exp(lvp)

    # 3. Micro-constants. Three transit transfers at ktr = 4/MTT and a
    #    terminal transfer into the central compartment at ka.
    ktr <- 4 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ka  * transit3
    d/dt(central)     <-  ka  * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # 5. Bioavailability applies to the dosing (GUT) compartment: F1 = FREL.
    f(depot) <- frel

    # 6. Observation, scaled with S5 = VC/1000 so that an amount in mg and a
    #    volume in L give a concentration in ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
