Tan_2009_artesunate <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model of oral artesunate (AS) and",
    "its active metabolite dihydroartemisinin (DHA) in 91 healthy Korean",
    "adult volunteers (Tan 2009). AS is described by a one-compartment",
    "first-order absorption / first-order elimination disposition; DHA by a",
    "two-compartment disposition (central + peripheral). AS is converted",
    "mole-for-mole to DHA as the only elimination pathway. Body weight",
    "linearly increases DHA apparent clearance (1.9 L/h per kg above the",
    "61.5 kg reference) and a high-fat / high-caloric meal at dosing reduces",
    "AS absorption rate Ka by 84%. Subjects pooled across four Phase I",
    "studies (single-dose ascending, drug-interaction with pyronaridine,",
    "food-effect, and three-day multiple-dose) at 2-5 mg/kg oral AS."
  )
  reference <- paste(
    "Tan B, Naik H, Jang IJ, Yu KS, Kirsch LE, Shin CS, Craft JC, Fleckenstein L.",
    "Population pharmacokinetics of artesunate and dihydroartemisinin following",
    "single- and multiple-dosing of oral artesunate in healthy subjects.",
    "Malaria Journal 2009;8:304. doi:10.1186/1475-2875-8-304"
  )
  vignette  <- "Tan_2009_artesunate"
  units     <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear-deviation effect on DHA apparent clearance, centred on the",
        "reference weight 61.5 kg (the combined-cohort median; Tan 2009",
        "Results p.7 / Table 1): cl_dha_typ = exp(lcl_dha) + e_wt_cl_dha *",
        "(WT - 61.5), with e_wt_cl_dha = 1.9 L/h per kg (Table 2).",
        "Time-fixed at baseline in the source analysis."
      ),
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-vs-fasted indicator at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "1 = subject received the AS dose 30 min after a high-fat /",
        "high-caloric test meal (food-effect study fed arm), 0 = fasted",
        "at dosing (the reference, used in the single-dose, drug-interaction,",
        "multiple-dose, and food-effect fasted arms). Multiplicative effect",
        "on the AS absorption rate Ka: ka_typ = exp(lka) * (1 + e_fed_ka *",
        "FED), with e_fed_ka = -0.84 (Tan 2009 Table 2; the high-fat meal",
        "reduces Ka by 84%, lengthening the absorption half-life from",
        "10.8 min to 67.5 min, per Discussion p.10). Per-record dose-state",
        "indicator: a multi-period subject can carry FED = 0 in one period",
        "and FED = 1 in another."
      ),
      source_name        = "FOOD"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 91L,
    n_studies       = 4L,
    n_profiles      = 118L,
    n_observations  = list(AS = 916L, DHA = 1352L),
    age_range       = "19-40 years (combined median 23)",
    weight_range    = "50.1-70 kg (combined median 61.5)",
    weight_median   = "61.5 kg (population reference for the WT effect on DHA CL/F)",
    sex_female_pct  = round(100 * (91 - 53) / 91, 1),
    race_ethnicity  = "Korean",
    disease_state   = "Healthy adult volunteers (no malaria)",
    dose_range      = "2-5 mg/kg oral artesunate, single dose or once daily for 3 days, with or without high-fat meal, with or without co-administered pyronaridine in 3:1 ratio",
    regions         = "Republic of Korea (Seoul National University Clinical Trial Center)",
    notes           = paste(
      "Demographics from Tan 2009 Table 1 (n=91 across the four Phase I",
      "sub-studies SP-C001-03: single-dose 28, drug-interaction with",
      "pyronaridine 19, food-effect 20, multiple-dose 24).",
      "9 subjects contributed two profiles in the drug-interaction study",
      "and 18 contributed two profiles in the food-effect study; the",
      "118 total profiles were treated as independent because the 21-day",
      "washout was much longer than the AS / DHA elimination half-lives",
      "(~12 min / ~40 min)."
    )
  )

  ini({
    # Structural parameters (Tan 2009 Table 2, "Estimate" column). The
    # paper fit AS and DHA simultaneously on natural-log nmol/L
    # concentrations with AS dose pre-converted to nmols (Methods p.4);
    # apparent volumes (L) and clearances (L/h) are unit-invariant
    # between mass and molar amount, so the published L and L/h values
    # apply directly when this model file tracks amounts in mg. The
    # mole-for-mole AS -> DHA conversion enters as a mass-rate factor
    # mw_dha / mw_ars at the metabolite-formation step (see model()).
    lka     <- log(3.85)
    label("Apparent first-order absorption rate constant for AS, Ka (1/h)")  # Tan 2009 Table 2: Ka = 3.85 1/h (%RSE 3.61) -- typical value at FED = 0
    lcl     <- log(1190)
    label("Apparent artesunate elimination clearance, CL/F (L/h)")  # Tan 2009 Table 2: CL/F = 1190 L/h (%RSE 4.20)
    lvc     <- log(1210)
    label("Apparent artesunate central volume of distribution, V2/F (L)")  # Tan 2009 Table 2: V2/F = 1210 L (%RSE 5.77)
    lcl_dha <- log(93.7)
    label("Apparent DHA elimination clearance from DHA central at WT = 61.5 kg, CLM/F (L/h)")  # Tan 2009 Table 2: CLM/F = 93.7 L/h (%RSE 3.30) -- typical value at the 61.5 kg reference weight
    lvc_dha <- log(97.1)
    label("Apparent DHA central volume of distribution, V3/F (L)")  # Tan 2009 Table 2: V3/F = 97.1 L (%RSE 4.85)
    lq_dha  <- log(5.74)
    label("Apparent DHA inter-compartmental clearance, Q/F (L/h)")  # Tan 2009 Table 2: Q/F = 5.74 L/h (%RSE 12.8)
    lvp_dha <- log(18.5)
    label("Apparent DHA peripheral volume of distribution, V4/F (L)")  # Tan 2009 Table 2: V4/F = 18.5 L (%RSE 10.6)

    # Covariate effects (Tan 2009 Results p.7).
    #
    # Food effect on AS absorption rate Ka. Paper text:
    #   "Ka = 3.85 * (1 + (-0.84) * FOOD)", with FOOD = 1 for fed and
    #   FOOD = 0 for fasted. The high-fat / high-caloric meal reduces Ka
    #   by 84% (Discussion p.10).
    # Encoded as a multiplicative deviation centred on FED = 0:
    #   ka_typ = exp(lka) * (1 + e_fed_ka * FED)
    # At FED = 1, multiplier = 1 - 0.84 = 0.16 (a 6.25-fold longer
    # absorption half-life, 10.8 min -> 67.5 min).
    e_fed_ka <- -0.84
    label("Effect of fed state on AS Ka (multiplicative; -0.84 = 84% reduction)")  # Tan 2009 Table 2: theta_FOOD-Ka = -0.84 (%RSE 2.32); Results p.7
    #
    # Body weight effect on DHA apparent clearance. Paper text:
    #   "CLM/F = 93.7 + 1.9 * (WT - 61.5)"
    # i.e. an absolute (linear) shift in L/h per kg of body weight
    # deviation from the 61.5 kg reference, not a power or multiplicative
    # effect. Encoded as:
    #   cl_dha_typ = exp(lcl_dha) + e_wt_cl_dha * (WT - 61.5)
    # so e_wt_cl_dha carries units of L/h per kg (= 1.9). For a 70 kg
    # subject, cl_dha_typ = 93.7 + 1.9 * 8.5 = 109.9 L/h.
    e_wt_cl_dha <- 1.9
    label("Linear effect of WT on DHA CL/F (L/h per kg above 61.5 kg reference)")  # Tan 2009 Table 2: theta_WT-CLM/F = 1.90 (%RSE 16.3); Results p.7

    # Inter-individual variability (Tan 2009 Table 2 "Estimate" column).
    # The source paper modelled IIV as log-normal multiplicative:
    #   P_i = P_pop * exp(eta_i)
    # and reported both the variance and the corresponding %CV
    # (approximated as sqrt(variance) per the paper's own footnote
    # "magnitude of IIV was expressed as coefficient of variation (%CV),
    # which was approximated by the square root of the variance
    # estimate"). The numeric values below are the variances on the eta
    # scale taken directly from Table 2; the comments show the published
    # %CV for cross-check (sqrt(variance) * 100). A diagonal covariance
    # matrix was used (Methods p.4: "the data did not support the
    # implementation of a full variance-covariance matrix"); Q/F and
    # V4/F had no estimable IIV (Results p.7: "Fixing the variance of
    # the random effects for Q/F and V4/F to zero had little influence
    # on the OFV"), so etalq_dha and etalvp_dha are absent.
    etalcl     ~ 0.131    # Tan 2009 Table 2: var(eta_CL_AS) = 0.131 (CV 36.2%, %RSE 17.8)
    etalvc     ~ 0.330    # Tan 2009 Table 2: var(eta_V2_AS) = 0.330 (CV 57.4%, %RSE 20.9)
    etalka     ~ 1.26     # Tan 2009 Table 2: var(eta_Ka) = 1.26 (CV 112%, %RSE 15.4)
    etalcl_dha ~ 0.0786   # Tan 2009 Table 2: var(eta_CL_DHA) = 0.0786 (CV 28.0%, %RSE 22.5)
    etalvc_dha ~ 0.0901   # Tan 2009 Table 2: var(eta_V3_DHA) = 0.0901 (CV 30.0%, %RSE 36.0)

    # Residual error. The source paper modelled natural-log plasma
    # concentrations with additive residual on the log scale (Methods
    # p.4: "C_ij = C_pred,ij + eps_ij" where C_ij is the log-transformed
    # observation per the preceding paragraph "concentrations were then
    # natural log-transformed before the analysis"). By the standing
    # convention in references/parameter-names.md, NONMEM additive-on-
    # log-scale residual maps to nlmixr2 proportional residual in linear
    # space, with propSd = SD on the log scale ~= CV in linear space to
    # first order. Table 2 reports RV as %CV (37.5% AS, 28.2% DHA), with
    # the corresponding log-scale variances 0.141 and 0.080.
    propSd     <- 0.375
    label("Proportional residual SD for artesunate plasma concentration")  # Tan 2009 Table 2: RV AS = 37.5% (variance 0.141, %RSE 9.73)
    propSd_dha <- 0.282
    label("Proportional residual SD for DHA plasma concentration")  # Tan 2009 Table 2: RV DHA = 28.2% (variance 0.080, %RSE 11.2)
  })

  model({
    # Molecular weights (g/mol) as reported in Tan 2009 Methods p.4
    # ("the molecular weights of AS and DHA are quite different (384.4
    # for AS and 284.9 for DHA)"). Used to apply the mole-for-mole AS
    # -> DHA conversion in mass-rate units, since this model file
    # tracks amounts in mg (oral mg/kg dosing) while the source paper
    # fit on molar amounts (nmol).
    mw_ars <- 384.4  # artesunate
    mw_dha <- 284.9  # dihydroartemisinin

    # Individual PK parameters with covariate effects.
    # AS Ka is multiplicatively reduced by 84% under a high-fat meal;
    # DHA CL/F is additively shifted by 1.9 L/h per kg of WT deviation
    # from the 61.5 kg reference.
    ka     <- exp(lka + etalka) * (1 + e_fed_ka * FED)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc     + etalvc)
    cl_dha <- (exp(lcl_dha) + e_wt_cl_dha * (WT - 61.5)) * exp(etalcl_dha)
    vc_dha <- exp(lvc_dha + etalvc_dha)
    q_dha  <- exp(lq_dha)
    vp_dha <- exp(lvp_dha)

    # Micro-constants. Under the paper's assumption of complete in-vivo
    # conversion of AS to DHA, the only elimination of AS is conversion
    # to DHA, so kel_AS = cl/vc. DHA is eliminated linearly from its
    # central compartment.
    kel     <- cl     / vc
    kel_dha <- cl_dha / vc_dha
    k12_dha <- q_dha  / vc_dha
    k21_dha <- q_dha  / vp_dha

    # ODE system. Oral dose targets depot. AS in central is converted
    # mole-for-mole to DHA central via the MW ratio (so a mg of AS
    # eliminated produces (mw_dha / mw_ars) = 0.7411 mg of DHA).
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot - kel * central
    d/dt(central_dha)      <-  kel * central * (mw_dha / mw_ars) -
                               kel_dha * central_dha -
                               k12_dha * central_dha +
                               k21_dha * peripheral1_dha
    d/dt(peripheral1_dha)  <-  k12_dha * central_dha -
                               k21_dha * peripheral1_dha

    # Plasma concentrations. With dose in mg and volumes in L, central
    # / vc is mg/L (numerically equal to ug/mL); multiply by 1000 for
    # ng/mL when comparing against the paper's reported concentration
    # ranges (LLOQ 1 ng/mL for both species per Methods p.3).
    Cc     <- central     / vc
    Cc_dha <- central_dha / vc_dha

    Cc     ~ prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
