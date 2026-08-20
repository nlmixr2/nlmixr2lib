Nagy_2017_obiltoxaximab <- function() {
  description <- paste(
    "Preclinical and clinical (NZW rabbit, cynomolgus macaque, human).",
    "Joint two-compartment population PK model for obiltoxaximab, the anti-protective-antigen",
    "chimeric IgG1(kappa) monoclonal antibody approved for inhalational anthrax under the US FDA",
    "Animal Rule. One file carries the three species-specific parameter sets used for",
    "animal-to-human dose translation (Nagy 2017 Supplementary Table S1), selected by the",
    "SPECIES_RABBIT / SPECIES_MACAQUE indicators with human as the reference species.",
    "Animals additionally have first-order intramuscular absorption (Ka, F1), and",
    "anthrax-infected subjects carry a parallel Michaelis-Menten elimination arm that",
    "approximates protective-antigen target-mediated drug disposition. For infected humans the",
    "macaque nonlinear-clearance component is carried over and allometrically scaled to human",
    "body size, exactly as the paper does to project infected-human exposures."
  )
  reference <- paste(
    "Nagy CF, Mondick J, Serbina N, Casey LS, Carpenter SE, French J, Guttendorf R.",
    "Animal-to-human dose translation of obiltoxaximab for treatment of inhalational anthrax",
    "under the US FDA animal rule.",
    "Clin Transl Sci. 2017;10(1):12-19. doi:10.1111/cts.12433.",
    "Parameter estimates are from Supplementary Table S1 (file CTS-10-12-s002);",
    "noncompartmental comparators are from Supplementary Tables S3 (CTS-10-12-s006)",
    "and S4 (CTS-10-12-s007).",
    "The companion survival model from the same paper is available as",
    "modellib('Nagy_2017_obiltoxaximab_survival').",
    sep = " "
  )
  vignette <- "Nagy_2017_obiltoxaximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor. The reference weight is species-specific:",
        "3.165 kg for NZW rabbits and 2.88 kg for cynomolgus macaques",
        "(Nagy 2017 Results, 'Animal and human population pharmacokinetic modeling'),",
        "and 75 kg for humans (Nagy 2017 Methods; the weight used for all published human",
        "simulations -- the paper does not state a separate human normalisation weight).",
        "Studied ranges: rabbits 2.9-4.0 kg and macaques 2.7-7.3 kg prior to challenge",
        "(Methods, 'Animal pharmacokinetic studies'); humans 50-125 kg (Discussion)."
      ),
      source_name        = "WT"
    ),
    SPECIES_RABBIT = list(
      description        = "New Zealand White rabbit indicator, 1 = NZW rabbit, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (human, the reference species)",
      notes              = paste(
        "Selects the NZW rabbit column of Nagy 2017 Supplementary Table S1.",
        "Mutually exclusive with SPECIES_MACAQUE; both 0 selects the human parameter set."
      ),
      source_name        = "SPECIES"
    ),
    SPECIES_MACAQUE = list(
      description        = "Cynomolgus macaque indicator, 1 = cynomolgus macaque, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (human, the reference species)",
      notes              = paste(
        "Selects the cynomolgus macaque column of Nagy 2017 Supplementary Table S1.",
        "Mutually exclusive with SPECIES_RABBIT; both 0 selects the human parameter set."
      ),
      source_name        = "SPECIES"
    ),
    DIS_ANTHRAX = list(
      description        = "Active inhalational anthrax infection indicator, 1 = infected, 0 = healthy / unexposed",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy, unexposed)",
      notes              = paste(
        "Switches on the parallel Michaelis-Menten elimination arm approximating",
        "protective-antigen target-mediated drug disposition. Nagy 2017 Results:",
        "'TMDD in infected NZW rabbits and cynomolgus macaques was approximated via parallel",
        "nonlinear elimination for infected animals only'. The human population PK model was fit",
        "to healthy volunteers only and therefore has no fitted nonlinear arm; Nagy 2017 Methods",
        "state that to predict infected-human exposures 'the nonlinear clearance model component",
        "from the cynomolgus macaque model was added to the human population PK model,",
        "allometrically scaled to human body size'. Setting DIS_ANTHRAX = 1 with both species",
        "indicators 0 reproduces that published infected-human projection.",
        "In the animal studies infection was aerosol challenge with a target 200 LD50 of",
        "Bacillus anthracis (Ames strain) spores."
      ),
      source_name        = "INFECTED"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "obiltoxaximab", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "obiltoxaximab", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "obiltoxaximab", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human + rabbit (New Zealand White) + cynomolgus macaque",
    n_subjects     = 758L,
    n_studies      = 10L,
    weight_range   = "rabbits 2.9-4.0 kg; macaques 2.7-7.3 kg; humans 50-125 kg",
    disease_state  = paste(
      "Pooled healthy and Bacillus anthracis (Ames) aerosol-challenged NZW rabbits and",
      "cynomolgus macaques, plus healthy human volunteers. Healthy and infected animal data",
      "were fit simultaneously within each species to describe disease effects on PK.",
      "No human anthrax patients were studied -- infected-human exposure is a model projection",
      "(this is the core of the FDA Animal Rule dose-translation argument)."
    ),
    dose_range     = paste(
      "Rabbits: 3, 10, 30 mg/kg i.v. and 10 mg/kg i.m. (healthy); 1, 4, 8, 16 mg/kg i.v.",
      "(infected). Macaques: 3, 10, 30 mg/kg i.v. and 10 mg/kg i.m. (healthy);",
      "4, 8, 16, 32 mg/kg i.v. (infected). Humans: 120-360 mg and 4, 8, 16 mg/kg i.v.",
      "infused over 90 min."
    ),
    notes          = paste(
      "n_studies = 10 per Nagy 2017 Methods ('PK data from 10 studies (two studies in rabbits,",
      "five studies in macaques, and three studies in healthy humans) were used for population",
      "PK model fitting'). The combined data set consisted of 791, 929, and 2,830 observations",
      "for rabbits, macaques, and humans, respectively. n_subjects = 758 is the sum of the",
      "enrolment counts in Table 1 for the ten model-fitting studies",
      "(rabbits 50 + 70; macaques 24 + 44 + 48 + 48 + 51; humans 45 + 108 + 280);",
      "Table 1 reports N per study rather than a pooled analysis-set size.",
      "Human demographics (age, sex, race) are not broken out in this paper;",
      "they are reported in the companion five-study safety/PK paper",
      "(Nagy 2016, Clin Ther 38:2083-2097)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Human parameter set (reference species) -- Nagy 2017 Supplementary
    # Table S1, "Humans" column. Typical values apply at the 75 kg weight
    # used for every published human simulation (Methods).
    # ------------------------------------------------------------------
    lcl <- log(0.233); label("Clearance, human (L/day)")                            # Table S1 humans CL = 0.233 (95% CI 0.224, 0.242)
    lvc <- log(3.21);  label("Central volume of distribution, human (L)")           # Table S1 humans Vc = 3.21 (95% CI 3.12, 3.29)
    lvp <- log(2.73);  label("Peripheral volume of distribution, human (L)")        # Table S1 humans Vp = 2.73 (95% CI 2.62, 2.84)
    lq  <- log(0.473); label("Intercompartmental clearance, human (L/day)")         # Table S1 humans Q  = 0.473 (95% CI 0.427, 0.516)

    # ------------------------------------------------------------------
    # NZW rabbit parameter set -- Table S1, "NZW Rabbits" column.
    # Typical values apply at the 3.165 kg reference weight.
    # ------------------------------------------------------------------
    lcl_rabbit <- log(0.0263); label("Clearance, NZW rabbit (L/day)")                     # Table S1 rabbit CL = 0.0263 (95% CI 0.0239, 0.0375)
    lvc_rabbit <- log(0.114);  label("Central volume of distribution, NZW rabbit (L)")    # Table S1 rabbit Vc = 0.114 (95% CI 0.104, 0.120)
    lvp_rabbit <- log(0.0744); label("Peripheral volume of distribution, NZW rabbit (L)") # Table S1 rabbit Vp = 0.0744 (95% CI 0.0598, 0.0970)
    lq_rabbit  <- log(0.119);  label("Intercompartmental clearance, NZW rabbit (L/day)")  # Table S1 rabbit Q  = 0.119 (95% CI 0.0830, 0.154)
    lka_rabbit <- log(0.961);  label("First-order i.m. absorption rate, NZW rabbit (1/day)") # Table S1 rabbit Ka = 0.961 (95% CI 0.270, 1.43)
    lfdepot_rabbit <- log(0.899); label("Absolute i.m. bioavailability, NZW rabbit (fraction)") # Table S1 rabbit F1 = 0.899 (95% CI 0.788, 1.43)

    # ------------------------------------------------------------------
    # Cynomolgus macaque parameter set -- Table S1, "Cynomolgus Macaques"
    # column. Typical values apply at the 2.88 kg reference weight.
    # ------------------------------------------------------------------
    lcl_macaque <- log(0.0191); label("Clearance, cynomolgus macaque (L/day)")                     # Table S1 macaque CL = 0.0191 (95% CI as printed 0.0162, 0.223; see vignette Errata)
    lvc_macaque <- log(0.134);  label("Central volume of distribution, cynomolgus macaque (L)")    # Table S1 macaque Vc = 0.134 (95% CI 0.127, 0.141)
    lvp_macaque <- log(0.123);  label("Peripheral volume of distribution, cynomolgus macaque (L)") # Table S1 macaque Vp = 0.123 (95% CI 0.109, 0.138)
    lq_macaque  <- log(0.0890); label("Intercompartmental clearance, cynomolgus macaque (L/day)")  # Table S1 macaque Q  = 0.0890 (95% CI 0.0785, 0.103)
    lka_macaque <- log(3.89);   label("First-order i.m. absorption rate, cynomolgus macaque (1/day)") # Table S1 macaque Ka = 3.89 (95% CI 2.96, 5.08)
    lfdepot_macaque <- log(0.895); label("Absolute i.m. bioavailability, cynomolgus macaque (fraction)") # Table S1 macaque F1 = 0.895 (95% CI 0.777, 1.02)

    # ------------------------------------------------------------------
    # Parallel Michaelis-Menten (protective-antigen TMDD surrogate) arm,
    # infected animals only. Vmax is encoded as an AMOUNT rate (mg/day).
    #
    # Nagy 2017 Table S1 labels the Vmax unit "ug/mL/day" (a concentration
    # rate). That unit label is not consistent with the paper's own
    # results; the value read as an amount rate (mg/day) is. Three
    # independent checks in the source agree on the amount-rate reading
    # and reject the concentration-rate reading -- see the vignette
    # "Assumptions and deviations" section for the worked numbers:
    #   (1) infected-macaque apparent CL at 16 mg/kg: 8.3 mL/day/kg
    #       (amount rate) vs 6.9 (concentration rate); Table S3 observed
    #       range is 7.8-9.0;
    #   (2) the stated 18% reduction in human AUC0-inf when TMDD is
    #       included (Results, Table 2): 18.1% (amount rate) vs 3.5%
    #       (concentration rate);
    #   (3) Table 2 infected-human AUC0-inf 4,070 ug*day/mL is reproduced
    #       only by the amount-rate reading.
    # Only the UNIT is reinterpreted; the numeric values are exactly as
    # published.
    # ------------------------------------------------------------------
    lvmax_rabbit  <- log(0.912); label("Maximum rate of PA-mediated elimination, NZW rabbit (mg/day)")          # Table S1 rabbit Vmax = 0.912 (95% CI 0.000506, 1.68)
    lkm_rabbit    <- log(10.4);  label("Obiltoxaximab concentration at half-maximal PA-mediated elimination, NZW rabbit (ug/mL)") # Table S1 rabbit Km = 10.4 (95% CI -0.266, 39.5)
    lvmax_macaque <- log(0.275); label("Maximum rate of PA-mediated elimination, cynomolgus macaque (mg/day)")  # Table S1 macaque Vmax = 0.275 (95% CI 0.0671, 0.903)
    lkm_macaque   <- log(3.21);  label("Obiltoxaximab concentration at half-maximal PA-mediated elimination, cynomolgus macaque (ug/mL)") # Table S1 macaque Km = 3.21 (95% CI 1.12, 35.6)

    # ------------------------------------------------------------------
    # Allometric exponents. Nagy 2017 states only that "volume and
    # clearance parameters were allometrically scaled, normalized to a
    # reference weight of 3.165 kg for NZW rabbits and 2.88 kg for
    # cynomolgus macaques" (Results) and does not report the exponents.
    # The standard values are assumed here.
    # ------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless; assumed standard value, not reported by Nagy 2017)")
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless; assumed standard value, not reported by Nagy 2017)")

    # ------------------------------------------------------------------
    # Nagy 2017 Supplementary Table S1 reports typical (fixed-effect)
    # values only. No IIV variances and no residual-error magnitudes are
    # published anywhere in the paper or its supplements, although the
    # visual predictive checks in Supplementary Figure S1 confirm both
    # were estimated. They are therefore fixed to zero rather than
    # invented; the model reproduces typical-value profiles only.
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; not published by Nagy 2017)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Species selection. Human is the reference species: both
    #    indicators 0 selects the human parameter set.
    # ------------------------------------------------------------------
    isRabbit  <- SPECIES_RABBIT
    isMacaque <- SPECIES_MACAQUE
    isHuman   <- 1 - SPECIES_RABBIT - SPECIES_MACAQUE

    # Species-specific allometric reference weight (kg): rabbit and
    # macaque values from Nagy 2017 Results; the human simulations were
    # all run at 75 kg (Methods).
    wtRef <- isRabbit * 3.165 + isMacaque * 2.88 + isHuman * 75

    # ------------------------------------------------------------------
    # 2. Species-specific typical values (Table S1), then allometric
    #    scaling to the subject's own body weight.
    # ------------------------------------------------------------------
    clRef <- isRabbit * exp(lcl_rabbit) + isMacaque * exp(lcl_macaque) + isHuman * exp(lcl)
    vcRef <- isRabbit * exp(lvc_rabbit) + isMacaque * exp(lvc_macaque) + isHuman * exp(lvc)
    vpRef <- isRabbit * exp(lvp_rabbit) + isMacaque * exp(lvp_macaque) + isHuman * exp(lvp)
    qRef  <- isRabbit * exp(lq_rabbit)  + isMacaque * exp(lq_macaque)  + isHuman * exp(lq)

    cl <- clRef * (WT / wtRef)^e_wt_cl
    vc <- vcRef * (WT / wtRef)^e_wt_vc
    vp <- vpRef * (WT / wtRef)^e_wt_vc
    q  <- qRef  * (WT / wtRef)^e_wt_cl

    # ------------------------------------------------------------------
    # 3. Intramuscular absorption. Only the two animal species were dosed
    #    i.m. (Table 1), so Table S1 has no human Ka / F1. The macaque
    #    values are carried over for humans purely so the depot route
    #    stays defined; i.m. dosing in humans is NOT supported by this
    #    source and should not be simulated. Human dosing in Nagy 2017 is
    #    i.v. only, which bypasses the depot entirely.
    # ------------------------------------------------------------------
    ka <- isRabbit * exp(lka_rabbit) + (1 - isRabbit) * exp(lka_macaque)
    fdepot <- isRabbit * exp(lfdepot_rabbit) + (1 - isRabbit) * exp(lfdepot_macaque)

    # ------------------------------------------------------------------
    # 4. Parallel Michaelis-Menten elimination approximating
    #    protective-antigen TMDD, present only when DIS_ANTHRAX = 1.
    #    Humans inherit the macaque component (Methods), so the Vmax
    #    allometric reference weight for humans is the MACAQUE reference
    #    weight (2.88 kg) rather than 75 kg -- the component is scaled
    #    from where it was estimated. Km is a concentration and is not
    #    scaled.
    # ------------------------------------------------------------------
    wtRefVmax <- isRabbit * 3.165 + (1 - isRabbit) * 2.88
    vmaxRef <- isRabbit * exp(lvmax_rabbit) + (1 - isRabbit) * exp(lvmax_macaque)
    km      <- isRabbit * exp(lkm_rabbit)   + (1 - isRabbit) * exp(lkm_macaque)

    vmax <- DIS_ANTHRAX * vmaxRef * (WT / wtRefVmax)^e_wt_cl

    # ------------------------------------------------------------------
    # 5. Two-compartment disposition with first-order depot absorption.
    #    central is in mg and vc in L, so Cc is in mg/L == ug/mL, matching
    #    both units$concentration and the Km units in Table S1. vmax is an
    #    amount rate in mg/day, so the Michaelis-Menten term is already an
    #    amount rate and needs no volume conversion.
    #
    #    The central concentration is written out IN FULL inside the
    #    Michaelis-Menten term rather than via an intermediate variable. In
    #    the nlmixr2 model-function form, an intermediate that references an
    #    ODE state and is assigned before the d/dt() line evaluates to 0
    #    inside the derivative, which silently drops the whole nonlinear term
    #    and makes infected subjects identical to healthy ones. Writing
    #    `central / vc` inline reproduces the plain-rxode2 solution exactly.
    #    Do not refactor this into a named intermediate.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          (q / vc) * central +
                          (q / vp) * peripheral1 -
                          vmax * (central / vc) / (km + central / vc)
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
