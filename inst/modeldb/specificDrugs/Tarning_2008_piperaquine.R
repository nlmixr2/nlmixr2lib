# Population pharmacokinetic model of oral piperaquine in Burmese and Karen
# adults and children with uncomplicated Plasmodium falciparum malaria
# treated on the Thai-Myanmar border with dihydroartemisinin-piperaquine
# (Artekin) given as three or four doses (Tarning 2008, Antimicrobial
# Agents and Chemotherapy 52(3):1052-1061; doi:10.1128/aac.00955-07).

Tarning_2008_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in Burmese and Karen adults",
    "and children with uncomplicated Plasmodium falciparum malaria",
    "(Tarning 2008). Two-compartment disposition with first-order",
    "absorption (no lag) and elimination from the central compartment.",
    "Body weight is the only retained covariate: a linear (1 + theta *",
    "(WT - 48)) effect on apparent oral clearance CL/F and on apparent",
    "central volume of distribution Vc/F, centred on the cohort median",
    "of 48 kg. The combined three-dose and four-dose Artekin regimens",
    "were pooled; no treatment-regimen effect was retained. Exponential",
    "IIV on all five disposition / absorption parameters. Residual error",
    "is proportional in linear concentration space (the source paper fit",
    "an additive error on natural-log-transformed concentrations).",
    sep = " "
  )
  reference <- paste(
    "Tarning J, Ashley EA, Lindegardh N, Stepniewska K, Phaiphun L,",
    "Day NPJ, McGready R, Ashton M, Nosten F, White NJ (2008).",
    "Population pharmacokinetics of piperaquine after two different",
    "treatment regimens with dihydroartemisinin-piperaquine in",
    "patients with Plasmodium falciparum malaria in Thailand.",
    "Antimicrobial Agents and Chemotherapy 52(3):1052-1061.",
    "doi:10.1128/aac.00955-07.",
    sep = " "
  )
  vignette <- "Tarning_2008_piperaquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear covariate effect centred on the cohort median 48 kg",
        "(Tarning 2008 Table 2). CL/F = 66.0 * (1 + 0.0262 * (WT - 48))",
        "and Vc/F = 8660 * (1 + 0.0273 * (WT - 48)). Methods (page 1054):",
        "'continuous covariates were evaluated by including them in the",
        "model as linear, allometric, or hyperbolic maximum effect",
        "functions centered on the median value' and Results (page 1054):",
        "'A linear relationship between body weight and clearance and",
        "body weight and central volume of distribution gave the best",
        "fit to the data (2.6% and 2.7% increase in oral clearance and",
        "central volume of distribution, respectively, per kg of body",
        "weight increase from median weight). Use of an allometric or a",
        "hyperbolic maximum effect model for body weight on clearance",
        "and/or intercompartment clearance did not converge with the",
        "FOCE method.' Time-fixed at admission in the source data set.",
        "Cohort weight range 12 - 74 kg (Table 1); the authors caution",
        "(page 1055) that 'the covariate function provided by these data",
        "should not be extrapolated beyond the studied population",
        "demographics, since for children below 10 kg of body weight",
        "parameter estimates will be unreasonable'.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 98L,
    n_studies       = 1L,
    n_dp4           = 50L,
    n_dp3           = 48L,
    n_children      = 11L,
    age_range       = "3 - 55 years (pooled, Table 1)",
    age_median      = "25 years in both DP4 and DP3 arms (Table 1)",
    weight_range    = "12 - 74 kg (pooled across DP4 14 - 74 kg and DP3 12 - 59 kg, Table 1)",
    weight_median   = "47 kg DP4 / 51 kg DP3; pooled cohort median 48 kg used in the WT-covariate centring (Tables 1 - 2)",
    height_range    = "92 - 170 cm (pooled, Table 1)",
    sex_female_pct  = 39.8,
    race_ethnicity  = "Burmese or Karen ethnicity (study sites of the Shoklo Malaria Research Unit on the Thai-Myanmar border)",
    disease_state   = paste(
      "Uncomplicated symptomatic Plasmodium falciparum malaria.",
      "Exclusion criteria: pregnancy, lactation, >= 4% parasitised red",
      "blood cells, age below 1 year or above 65 years, signs or",
      "symptoms of severe malaria, or treatment with mefloquine within",
      "the previous 60 days (Methods, page 1053)."
    ),
    dose_range      = paste(
      "Dihydroartemisinin-piperaquine fixed-dose combination (Artekin;",
      "Holleykin Pharmaceutical Co. Ltd., Guangzhou, China; each tablet",
      "contains 320 mg piperaquine phosphate + 40 mg dihydroartemisinin)",
      "titrated to a total dose of 55 mg/kg piperaquine phosphate",
      "(equivalent to ~31 mg/kg piperaquine base; Table 1) rounded up to",
      "the nearest half tablet. Randomly assigned to either the standard",
      "four-dose regimen DP4 (split equally at 0, 8, 24, 48 h) or the",
      "once-daily three-dose regimen DP3 (split equally at 0, 24, 48 h)."
    ),
    regions         = "Thai-Myanmar border (Shoklo Malaria Research Unit, Mae Sot, Thailand)",
    notes           = paste(
      "Demographics from Tarning 2008 Table 1. 480 venous plasma samples",
      "drawn from 98 patients over 63 days post-treatment; 469 had",
      "measurable piperaquine concentrations and 11 samples below the",
      "lower limit of quantification (5 ng/mL or 2.5 ng/mL depending on",
      "the plasma aliquot volume) were coded as missing. Modelled in",
      "NONMEM v V level 1.1 with FOCE estimation on natural-log-",
      "transformed concentrations; visual predictive check used 500",
      "simulations per sampling time (Methods, page 1054)."
    )
  )

  ini({
    # Structural parameters from Tarning 2008 Table 2 ("Population estimate
    # (% RSE)" column). Apparent values (relative to F) reported on the
    # linear scale; log() applied here for the nlmixr2 internal log scale.
    # Typical values for CL/F and Vc/F correspond to the cohort median
    # weight of 48 kg (the centring value used by the WT covariate model).

    lcl <- log(66.0)
    label("Apparent piperaquine elimination clearance CL/F at WT = 48 kg (L/h)")
    # Tarning 2008 Table 2: CL/F = 66.0 L/h (% RSE 6.9) at WT = 48 kg

    lvc <- log(8660)
    label("Apparent central volume of distribution Vc/F at WT = 48 kg (L)")
    # Tarning 2008 Table 2: Vc/F = 8,660 L (% RSE 14) at WT = 48 kg

    lq <- log(131)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Tarning 2008 Table 2: Q/F = 131 L/h (% RSE 13)

    lvp <- log(24000)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Tarning 2008 Table 2: Vp/F = 24,000 L (% RSE 13)

    lka <- log(0.717)
    label("First-order absorption rate constant ka (1/h)")
    # Tarning 2008 Table 2: ka = 0.717 h^-1 (% RSE 25)

    # Linear WT covariate effects centred on the cohort median 48 kg.
    # Tarning 2008 Table 2 gives the covariate parameter as the per-kg
    # fractional increase from median weight in the multiplicative
    # bracket {1 + theta * (WT - 48)}, with the bracket applied to both
    # CL/F and Vc/F.
    e_wt_cl <- 0.0262
    label("Linear WT covariate effect on CL/F (per kg above 48 kg)")
    # Tarning 2008 Table 2: CL/F covariate = 0.0262 (% RSE 2.9) per kg

    e_wt_vc <- 0.0273
    label("Linear WT covariate effect on Vc/F (per kg above 48 kg)")
    # Tarning 2008 Table 2: Vc/F covariate = 0.0273 (% RSE 11) per kg

    # IIV. Tarning 2008 Table 2 reports inter-individual variability as
    # "% CV for interindividual variability (% RSE)" alongside each
    # parameter; the variability is modelled exponentially (Methods,
    # page 1054: "Interindividual random variability in all parameters
    # was modeled exponentially as illustrated for clearance:
    # (CL/F)_i = TV(CL/F) * exp(eta_i,CL/F) ..."). Following the
    # log-normal convention CV = sqrt(exp(omega^2) - 1), the internal
    # log-scale variance omega^2 = log(CV^2 + 1).
    #
    #   CL/F  IIV 42%   -> omega^2 = log(0.42^2 + 1) = 0.162519
    #   Vc/F  IIV 101%  -> omega^2 = log(1.01^2 + 1) = 0.703107
    #   Q/F   IIV 85%   -> omega^2 = log(0.85^2 + 1) = 0.543763
    #   Vp/F  IIV 50%   -> omega^2 = log(0.50^2 + 1) = 0.223144
    #   ka    IIV 168%  -> omega^2 = log(1.68^2 + 1) = 1.340808
    etalcl ~ 0.162519
    # Tarning 2008 Table 2: IIV on CL/F = 42% CV (% RSE 44)

    etalvc ~ 0.703107
    # Tarning 2008 Table 2: IIV on Vc/F = 101% CV (% RSE 17)

    etalq ~ 0.543763
    # Tarning 2008 Table 2: IIV on Q/F = 85% CV (% RSE 18)

    etalvp ~ 0.223144
    # Tarning 2008 Table 2: IIV on Vp/F = 50% CV (% RSE 76)

    etalka ~ 1.340808
    # Tarning 2008 Table 2: IIV on ka = 168% CV (% RSE 38)

    # Residual error. Methods (page 1054): piperaquine plasma
    # concentrations "were transformed into their natural logarithms"
    # before modelling; the Discussion (page 1055) clarifies that
    # "the effective random residual error model should be considered
    # multiplicative, since the modeled data was log transformed."
    # Table 2 reports sigma = 31.4% CV, i.e. the proportional CV in
    # linear concentration space; mapped directly to nlmixr2 propSd.
    propSd <- 0.314
    label("Proportional residual SD for piperaquine plasma concentration (fraction)")
    # Tarning 2008 Table 2: sigma = 31.4% CV (% RSE 29)
  })

  model({
    # Individual PK parameters. Exponential IIV around the typical
    # values; WT covariate effects applied multiplicatively to CL/F and
    # Vc/F in the linear (1 + theta * (WT - 48)) form (Tarning 2008
    # Table 2). Q/F, Vp/F, and ka are not covariate-adjusted.
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 48))
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 48))
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order absorption from a depot compartment into a
    # two-compartment disposition model with elimination from central
    # (Tarning 2008 Methods, page 1054).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Piperaquine plasma concentration. Dose in mg of piperaquine base
    # (Tarning 2008 Table 1 reports the per-patient total piperaquine
    # dose in mg/kg base alongside the per-tablet piperaquine phosphate
    # content; the Table 2 CL/F = 66 L/h and the comparison column in
    # Table 3 give 1.4 L/h/kg, consistent with the base-mass-equivalent
    # parameters used in the Hung 2004 and other compared studies in
    # Table 3). Vc in L gives central / vc in mg/L; the 1000 factor
    # converts to ng/mL, matching the LOQ and assay-quality-control
    # concentrations reported in Methods (page 1053).
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale
    # (see ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
