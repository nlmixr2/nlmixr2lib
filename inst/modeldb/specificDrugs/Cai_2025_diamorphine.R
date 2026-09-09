Cai_2025_diamorphine <- function() {
  description <- paste(
    "Integrated four-compartment parent-plus-two-metabolite population PK",
    "model for diamorphine (heroin) and its sequential deacetylation products",
    "6-monoacetylmorphine (6-MAM) and morphine, after intramuscular or",
    "intranasal diamorphine in adult male heroin users (Cai 2025).",
    "Diamorphine reaches a one-compartment central pool through a common",
    "first-order absorption rate constant shared by both routes; the",
    "intramuscular route is the reference (F = 1) and the intranasal route",
    "carries an estimated relative bioavailability of about 52% on a",
    "logit-normal scale. Diamorphine and 6-MAM have no elimination other than",
    "sequential conversion (the paper discards their renal elimination and",
    "treats both as fully converted to morphine), so each is described by a",
    "single first-order rate constant; the 6-MAM central volume is set equal",
    "to the diamorphine central volume because it is not identifiable.",
    "Morphine follows a two-compartment disposition with its own central",
    "volume and clearance. All states are molar amounts, so the sequential",
    "1:1 deacetylation transfers amounts directly. Every parameter is",
    "standardised to 70 kg by theory-based allometry with fixed exponents",
    "(0.75 for clearance, 1 for volumes, -0.25 for first-order rate",
    "constants), and morphine clearance additionally carries a sigmoid",
    "postmenstrual-age maturation function fixed from Holford 2012",
    "(TM50 = 58.1 weeks, Hill = 3.58) that allows extrapolation from the",
    "adult fit to children."
  )
  reference <- paste(
    "Cai L, Zhai J, Ji B, Han F, Niu T, Wang L, Wang J.",
    "Intranasal diamorphine population pharmacokinetics modeling and",
    "simulation in pediatric breakthrough pain.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(3):435-447.",
    "doi:10.1002/psp4.13186.",
    sep = " "
  )
  vignette <- "Cai_2025_diamorphine"

  # Both administration routes are dosing targets: `depot` is intramuscular
  # (the reference, F = 1) and `depot2` is intranasal. Declared explicitly
  # because the automatic detection only recognises `depot` and `central`.
  dosing <- c("depot", "depot2")

  # The paper converted every measured concentration from ng/mL to nM
  # (Methods, "Dataset preparation") so that the sequential deacetylation of
  # diamorphine (369.4 g/mol) to 6-MAM (327.4 g/mol) to morphine
  # (285.34 g/mol) transfers molar amounts 1:1 without a molecular-weight
  # ratio on any transfer arrow. States are therefore nmol and observations
  # nmol/L; a mass dose must be converted before use (see the vignette).
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482: what each ODE state holds. Every state is a molar amount of
  # the named analyte, in nmol. The two depots are the two administration
  # routes the paper modelled as separate `depot()` macros both targeting the
  # diamorphine central compartment (Supplementary Material S1, [LONGITUDINAL]
  # PK block).
  compartmentData <- list(
    depot = list(
      analyte = "diamorphine", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    depot2 = list(
      analyte = "diamorphine", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "diamorphine", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_6mam = list(
      analyte = "6-monoacetylmorphine", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_morphine = list(
      analyte = "morphine", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_morphine = list(
      analyte = "morphine", units = "nmol",
      specimen = "tissue", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters every structural parameter a priori",
        "through theory-based allometry standardised to 70 kg",
        "(Supplementary Material S2, Eqs. 1-3): clearance scales as",
        "(WT/70)^0.75, volumes as (WT/70)^1 and first-order rate constants",
        "as (WT/70)^-0.25. The exponents are fixed, not estimated -- the",
        "supplement states the weight influence was 'fixed via empirical",
        "allometric scaling functions', and the main text (Results,",
        "'Statistical modeling results') reports that estimating the",
        "exponents worsened the fit. The rate-constant exponent is exactly",
        "the clearance exponent minus the volume exponent because a",
        "first-order rate constant is a clearance divided by a volume."
      ),
      source_name        = "WT (Supplementary Material S1 [COVARIATE] block)"
    ),
    AGE = list(
      description        = "Postnatal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Used only to derive postmenstrual age for",
        "the morphine-clearance maturation function:",
        "PMA (weeks) = 40 + AGE * 52 (Supplementary Material S2, Eq. 4).",
        "PMA then enters morphine clearance as the sigmoid maturation",
        "multiplier 1 / (1 + (TM50 / PMA)^Hill) with TM50 = 58.1 weeks and",
        "Hill = 3.58 (Eq. 5), both fixed from Holford 2012. The model was",
        "fitted only to adults, in whom the multiplier is numerically 1; the",
        "maturation term exists so the adult fit can be extrapolated to",
        "children, which is the paper's purpose. Because PMA is derived",
        "inside model() the data need carry only AGE, matching the paper's",
        "own covariate input set {AGE, WT}."
      ),
      source_name        = "AGE (Supplementary Material S1 [COVARIATE] block)"
    )
  )

  # Screened in the stepwise covariate modelling and NOT retained in the
  # final model (Results, "Statistical modeling results": "The automatic
  # stepwise covariates modeling did not identify other significant
  # covariates"). Recorded here so the paper's covariate screen is preserved
  # without declaring entries that model() never references. Note that the
  # administration route still enters the model structurally, as the relative
  # bioavailability of the intranasal depot; what the screen rejected was an
  # additional route effect on the absorption rate constant (Supplementary
  # Material S2: "we only assessed the route effect on diamorphine
  # absorption").
  covariatesDataExcluded <- list(
    DOSE_HIGH = list(
      description = "Diamorphine dose level indicator (0 = 6 mg, 1 = 12 mg)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "CAT2 in Supplementary Material S2, Eqs. 9-10 and Table S2. Screened",
        "on the disposition parameters and rejected, which the paper reads as",
        "confirming linear kinetics over the twofold dose range studied",
        "(Results: 'diamorphine absorption and disposition were not",
        "dosage-dependent')."
      )
    ),
    ROUTE_IN = list(
      description = "Administration route indicator (0 = intramuscular, 1 = intranasal)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "CAT1 in Supplementary Material S2, Eqs. 9-10. Screened as an",
        "exponential effect on the absorption rate constant and rejected:",
        "'When we separately estimated the absorption rate constant (Ka) for",
        "the two administration routes, there were minor differences, and",
        "using a common Ka improved the model fit' (Results). The route is",
        "still represented in the model, but structurally -- by which depot",
        "the dose enters -- rather than as a covariate."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10,
    n_studies      = 2,
    age_range      = "23-41 years",
    age_median     = "28.25 years (IM 6 mg) / 31.38 years (IN 6 and 12 mg)",
    weight_range   = "60.4-81.4 kg",
    weight_median  = "73.19 kg (IM 6 mg) / 72.16 kg (IN 6 and 12 mg)",
    sex_female_pct = 0,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Regular heroin users studied after at least three consecutive days of",
      "negative opioid tests (abstinence), i.e. opioid-free at dosing. Not a",
      "pain population: the paper fits adults and then extrapolates to",
      "children with breakthrough pain."
    ),
    dose_range     = paste(
      "Single doses of diamorphine hydrochloride 6 mg intramuscular",
      "(reference) and 6 or 12 mg intranasal, in a double-blind",
      "double-dummy crossover with a 1-week washout. 6 mg is roughly",
      "0.07-0.09 mg/kg in this cohort."
    ),
    regions        = "United States (NIDA Addiction Research Center, Baltimore) and Germany",
    notes          = paste(
      "Development dataset pooled from two published NIDA studies sharing",
      "dosing regimen, study medication, assay and sampling times",
      "(Cone et al. and Skopp et al.; paper Table 1). 385 plasma",
      "concentrations of diamorphine, 6-MAM and morphine by GC-MS.",
      "Demographics from Table S2. Between-occasion variability could not be",
      "estimated (unsuccessful minimisation or highly imprecise estimates),",
      "so each of the 28 treatment sessions across the 10 volunteers was",
      "treated as one 'modeled subject' and all random effects are",
      "between-subject variability on that basis (Results, 'Statistical",
      "modeling results'). Two further studies were used only for external",
      "verification and contributed no parameter estimates: Girardin et al.",
      "(8 adults, IM 67-202 mg) and Kidd et al. (12 children aged 4-13 y,",
      "IN 0.1 mg/kg). Estimation used Monolix 2021R2 (SAEM) with Bayesian",
      "priors on F, V1, V3 and CL; all parameter values encoded here are the",
      "final Table 2 estimates, not the Table S1 initial estimates.",
      "Females were excluded from the analysis (paper Discussion,",
      "limitations)."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Structural parameters. All typical values are standardised to a
    # 70 kg adult; Table 2 column "Estimation (RSE%)".
    #
    # Naming note: diamorphine and 6-MAM have no elimination pathway other
    # than sequential conversion, so the paper's metabolic transfer rate
    # constants K12 (diamorphine -> 6-MAM) and K23 (6-MAM -> morphine) ARE
    # the total elimination rate constants of those two species and take the
    # canonical `kel` names. The paper's K3p / Kp3 are the morphine
    # central <-> peripheral micro-constants.
    # ----------------------------------------------------------------
    lka <- log(3.04)
    label("Absorption rate constant, shared by the IM and IN routes (1/h)")  # Table 2: Ka = 3.04 (RSE 12.9%); a common Ka fitted better than route-specific values

    logitfdepot <- log(0.519 / (1 - 0.519))
    label("Logit of intranasal relative bioavailability versus intramuscular (fraction)")  # Table 2: F% = 0.519 (RSE 13.6%); logit(0.519) = 0.07603; IM is the reference route with F = 1

    lvc <- log(8.21)
    label("Diamorphine central volume at 70 kg, also used for 6-MAM (L)")  # Table 2: V1 = 8.21 L/70 kg (RSE 28.9%); V2 = V1 because the 6-MAM volume is unidentifiable

    lkel <- log(103)
    label("Diamorphine elimination rate constant, all of which forms 6-MAM (1/h)")  # Table 2: K12 = 103 /h (RSE 23.5%)

    lkel_6mam <- log(106)
    label("6-MAM elimination rate constant, all of which forms morphine (1/h)")  # Table 2: K23 = 106 /h (RSE 23.3%)

    lvc_morphine <- log(32.5)
    label("Morphine central volume at 70 kg (L)")  # Table 2: V3 = 32.5 L/70 kg (RSE 13.5%)

    lcl_morphine <- log(132)
    label("Morphine clearance at 70 kg in a mature subject (L/h)")  # Table 2: CL = 132 L/h/70 kg (RSE 19.0%)

    lk12_morphine <- log(24.2)
    label("Morphine central-to-peripheral rate constant (1/h)")  # Table 2: K3p = 24.2 /h (RSE 14.1%)

    lk21_morphine <- log(2.69)
    label("Morphine peripheral-to-central rate constant (1/h)")  # Table 2: Kp3 = 2.69 /h (RSE 14.3%)

    # ----------------------------------------------------------------
    # Allometric exponents. Supplementary Material S2 Eqs. 1-3 fixes these
    # at the theory-based values; the main text records that estimating them
    # degraded parameter precision. Named with the shared-exponent form
    # because a single exponent is applied across every clearance-type and
    # every volume-type parameter in the model.
    # ----------------------------------------------------------------
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on morphine clearance (unitless)")  # Supplementary Material S2 Eq. 1

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on the diamorphine and morphine volumes (unitless)")  # Supplementary Material S2 Eq. 3

    # ----------------------------------------------------------------
    # Morphine clearance maturation. Fixed from Holford 2012 (the paper's
    # supplement reference 2); reported as constants with no uncertainty.
    # The model has exactly one clearance, so these need no analyte suffix.
    # ----------------------------------------------------------------
    ltm50_cl <- fixed(log(58.1))
    label("Postmenstrual age at which morphine clearance reaches 50% of the mature value, literature value (weeks)")  # Supplementary Material S2 Eq. 5: TM50 = 58.1

    e_age_cl_hill <- fixed(3.58)
    label("Hill coefficient of the morphine clearance maturation function, literature value (unitless)")  # Supplementary Material S2 Eq. 5: Hill_CL = 3.58

    # ----------------------------------------------------------------
    # Between-subject variability. Table 2 column "BSV (RSE%)". Monolix
    # reports these as standard deviations on the transformed scale
    # (Supplementary Material S1 declares each parameter with `sd=omega_*`,
    # and Table S1 notes the default OMEGA initial value of 1 was used),
    # so each variance below is the squared table entry.
    # ----------------------------------------------------------------
    etalka + etalvc ~ c(0.609^2, 0.248601, 0.478^2)  # Table 2: BSV Ka = 0.609, BSV V1 = 0.478, corr_V1_Ka = 0.854; covariance = 0.854 * 0.609 * 0.478

    etalogitfdepot  ~ 0.568^2  # Table 2: BSV F% = 0.568, on the logit scale
    etalkel         ~ 0.400^2  # Table 2: BSV K12 = 0.400
    etalkel_6mam    ~ 0.295^2  # Table 2: BSV K23 = 0.295
    etalvc_morphine ~ 0.279^2  # Table 2: BSV V3 = 0.279
    etalcl_morphine ~ 0.297^2  # Table 2: BSV CL = 0.297
    etalk12_morphine ~ 0.125^2  # Table 2: BSV K3p = 0.125
    etalk21_morphine ~ 0.385^2  # Table 2: BSV Kp3 = 0.385

    # ----------------------------------------------------------------
    # Residual unexplained variability. Proportional error on each of the
    # three analytes (Results: "Proportional error models were employed to
    # account for RUV"; Supplementary Material S1 declares
    # errorModel=proportional(b1..b3) on C1, C2 and C3).
    # ----------------------------------------------------------------
    propSd <- 0.430
    label("Proportional residual error for diamorphine (fraction)")  # Table 2: RUV1 = 0.430 (RSE 12.0%)

    propSd_6mam <- 0.215
    label("Proportional residual error for 6-MAM (fraction)")  # Table 2: RUV2 = 0.215 (RSE 8.34%)

    propSd_morphine <- 0.236
    label("Proportional residual error for morphine (fraction)")  # Table 2: RUV3 = 0.236 (RSE 6.74%)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Derived covariate terms.
    # ----------------------------------------------------------------
    wt_ref <- 70  # allometric reference weight (kg), Supplementary Material S2 Eqs. 1-3

    allom_cl <- (WT / wt_ref)^e_wt_cl_q   # Eq. 1
    allom_v  <- (WT / wt_ref)^e_wt_vc_vp  # Eq. 3

    # Eq. 2 gives the first-order rate-constant exponent as -0.25. Under
    # theory-based allometry that is exactly the clearance exponent minus the
    # volume exponent, because a rate constant is a clearance divided by a
    # volume: 0.75 - 1 = -0.25. It is derived from the two fixed exponents
    # rather than declared a third time, which keeps the three scalings
    # mutually consistent.
    allo_k <- (WT / wt_ref)^(e_wt_cl_q - e_wt_vc_vp)

    # Postmenstrual age in weeks from postnatal age in years, Eq. 4, and the
    # sigmoid maturation multiplier on morphine clearance, Eq. 5. Eq. 5 is
    # printed as 1 / (1 + (PMA/TM50)^-Hill), which is the same function as
    # the supplement's Mlxtran code 1 / (1 + (TM50/PMA)^Hill); the latter
    # form is used here. In an adult the multiplier is numerically 1.
    pma           <- 40 + AGE * 52
    tm50_cl       <- exp(ltm50_cl)
    maturation_cl <- 1 / (1 + (tm50_cl / pma)^e_age_cl_hill)

    # ----------------------------------------------------------------
    # 2. Individual parameters. Every parameter carries BSV (Table 2) and
    #    the allometric term appropriate to its dimension.
    # ----------------------------------------------------------------
    ka     <- exp(lka + etalka) * allo_k
    fdepot <- expit(logitfdepot + etalogitfdepot)
    vc     <- exp(lvc + etalvc) * allom_v

    kel      <- exp(lkel + etalkel) * allo_k
    kel_6mam <- exp(lkel_6mam + etalkel_6mam) * allo_k

    vc_morphine  <- exp(lvc_morphine + etalvc_morphine) * allom_v
    cl_morphine  <- exp(lcl_morphine + etalcl_morphine) * allom_cl * maturation_cl
    k12_morphine <- exp(lk12_morphine + etalk12_morphine) * allo_k
    k21_morphine <- exp(lk21_morphine + etalk21_morphine) * allo_k

    # ----------------------------------------------------------------
    # 3. Micro-constants.
    # ----------------------------------------------------------------
    kel_morphine <- cl_morphine / vc_morphine  # k30 = CL/V3 in Supplementary Material S1

    # ----------------------------------------------------------------
    # 4. ODE system, Supplementary Material S1 and paper Figure 1b.
    #    `depot` is the intramuscular route, modelled as the reference with
    #    F = 1; `depot2` is the intranasal route and carries the estimated
    #    relative bioavailability. Both absorb into the diamorphine central
    #    compartment with the same rate constant.
    #
    #    Every state is a molar amount, so the sequential deacetylation
    #    diamorphine -> 6-MAM -> morphine transfers amounts 1:1 with no
    #    molecular-weight ratio. Diamorphine and 6-MAM have no elimination
    #    other than that conversion: the paper discards their renal
    #    elimination and treats both as fully converted to morphine
    #    (Discussion).
    # ----------------------------------------------------------------
    d/dt(depot)  <- -ka * depot
    d/dt(depot2) <- -ka * depot2

    d/dt(central) <- ka * depot + ka * depot2 - kel * central

    d/dt(central_6mam) <- kel * central - kel_6mam * central_6mam

    d/dt(central_morphine) <-
      kel_6mam * central_6mam +
      k21_morphine * peripheral1_morphine -
      k12_morphine * central_morphine -
      kel_morphine * central_morphine

    d/dt(peripheral1_morphine) <-
      k12_morphine * central_morphine - k21_morphine * peripheral1_morphine

    # ----------------------------------------------------------------
    # 5. Bioavailability. The intramuscular depot is the reference group and
    #    is left at F = 1 (Methods: "we first modeled IMD data as the
    #    reference group and thus bioavailability, F%, is regarded as 100%").
    # ----------------------------------------------------------------
    f(depot2) <- fdepot

    # ----------------------------------------------------------------
    # 6. Observations, in nmol/L. The 6-MAM concentration uses the
    #    diamorphine volume because V2 = V1 (Supplementary Material S1:
    #    C2 = A2/V1).
    # ----------------------------------------------------------------
    Cc          <- central          / vc
    Cc_6mam     <- central_6mam     / vc
    Cc_morphine <- central_morphine / vc_morphine

    Cc          ~ prop(propSd)
    Cc_6mam     ~ prop(propSd_6mam)
    Cc_morphine ~ prop(propSd_morphine)
  })
}
