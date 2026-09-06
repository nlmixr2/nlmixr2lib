FernandezTeruel_2024_capivasertib <- function() {
  description <- "Three-compartment population PK model for capivasertib (oral pan-AKT inhibitor) with parallel first-order and zero-order absorption, absorption lag time, and sigmoidal time-dependent auto-inhibition of apparent clearance, in patients with advanced or metastatic solid tumours (Fernandez-Teruel 2024)"
  reference <- "Fernandez-Teruel C, Cullberg M, Eberlein C, Barry ST, Zhou D. Population pharmacokinetics of capivasertib in patients with advanced or metastatic solid tumours. Clin Pharmacokinet. 2024;63(9):1191-1201. doi:10.1007/s40262-024-01407-x"
  vignette <- "FernandezTeruel_2024_capivasertib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "capivasertib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BBW (baseline body weight). Reference 67 kg, the pooled",
        "median of Fernandez-Teruel 2024 Table 2. Enters twice, in two different",
        "functional forms taken verbatim from the paper's covariate equations:",
        "a LINEAR-deviation effect on the initial apparent clearance,",
        "CL0/F * (1 + (WT - 67) * e_wt_cl), and a POWER effect on the",
        "logit-scale first-order absorbed fraction,",
        "logitffo * (WT / 67)^e_wt_logitffo. Baseline-only (time-fixed);",
        "the analysis carried a single per-patient weight."
      ),
      source_name        = "BBW"
    ),
    FASTED_STRICT = list(
      description        = "Overnight-fasted dosing indicator (1 = dose taken after an overnight fast; 0 = semi-fasted or fed)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (semi-fasted or fed)",
      notes              = paste(
        "Source column FASTED. The paper's three prandial protocols collapse to",
        "a strict-vs-relaxed contrast, which is exactly this canonical's axis:",
        "an overnight fast is FASTED_STRICT = 1, while the semi-fasted protocol",
        "(no food intake from 2 h before to 1 h after the dose) and the fed state",
        "are pooled at FASTED_STRICT = 0. Same orientation as the source column,",
        "so no value transformation is needed. Gates the absorption lag time:",
        "there is no lag under an overnight fast (Fernandez-Teruel 2024",
        "Sect. 3.2), so the whole ALAG1 term is multiplied by",
        "(1 - FASTED_STRICT). Per dose record, not per subject: the OAK study",
        "dosed the same patients under both prandial states in a fixed-sequence",
        "crossover."
      ),
      source_name        = "FASTED"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (1 = capsule; 0 = tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = paste(
        "Source column CAP. The comparator is the capivasertib tablet, and both",
        "arms carry F = 1 in the apparent (/F) parameterisation: unlike most",
        "members of this family the formulation does NOT act on bioavailability",
        "or on Ka. Fernandez-Teruel 2024 Sect. 3.2 states 'no other differences",
        "were detected between tablets and capsules', and the only retained",
        "effect is on the absorption lag time (0.46 h for the capsule vs 0.212 h",
        "for the tablet, both under non-overnight-fasted conditions). Same",
        "narrow role as in Xu_2025_aficamten.R. Per dose record: the OAK study",
        "compared the two formulations within patient."
      ),
      source_name        = "CAP"
    ),
    DOSE_CAPIVASERTIB_MG = list(
      description        = "Planned capivasertib dose level per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column DOSE. The planned twice-daily dose level, 80-800 mg",
        "(Fernandez-Teruel 2024 Table 2 'Capivasertib planned dose'), not the",
        "daily total. Enters the maximal auto-inhibition of CL/F as a",
        "linear-deviation term centred on 480 mg,",
        "(1 + (DOSE_CAPIVASERTIB_MG - 480) * e_dose_cl_time_max), applied to the",
        "LOG-scale magnitude. Because e_dose_cl_time_max is negative and the",
        "log-scale magnitude is itself negative, the extent of auto-inhibition",
        "GROWS with dose: 21% at the 480 mg reference, 53% at 800 mg. Carried as",
        "a drug-specific column rather than a bare `DOSE` because rxode2's",
        "etTrans consumes a column literally named DOSE and never exposes it to",
        "model()."
      ),
      source_name        = "DOSE"
    ),
    CONMED_PACLITAXEL = list(
      description        = "Concomitant paclitaxel indicator (1 = capivasertib given with paclitaxel)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capivasertib monotherapy, or with fulvestrant)",
      notes              = paste(
        "Source column PACL. Paclitaxel was given only in the BEECH study",
        "(90 of 441 patients, 20.4%); the reference level pools capivasertib",
        "monotherapy and capivasertib + fulvestrant, since fulvestrant was not a",
        "significant covariate. Acts on the log-scale magnitude of the",
        "auto-inhibition and nearly abolishes it (a 2.15-fold multiplier on a",
        "negative log-scale value), giving about 20% higher steady-state CL/F",
        "with paclitaxel. Fernandez-Teruel 2024 Discussion cautions that",
        "paclitaxel is confounded with the BEECH study, whose sparse post-144 h",
        "sampling limits the ability to resolve the time-dependent PK."
      ),
      source_name        = "PACL"
    )
  )

  # Covariates that Fernandez-Teruel 2024 screened in the PsN stepwise
  # covariate-model build (forward p = 0.005, backward p = 0.001) and did NOT
  # retain in the final model. Documentation only -- the paper reports no point
  # estimate for any of them, so none can be encoded, and none is referenced in
  # model(). Sources: Fernandez-Teruel 2024 Sect. 2.4.3 (covariate list),
  # Sect. 3.3 / Abstract Results and Sect. 5 Conclusions (the not-significant
  # verdict).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened; not significant. Median 56 years (range 27-87), Table 2."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant. 79.6% female, Table 2."
    ),
    CRCL = list(
      description = "Creatinine clearance", units = "mL/min", type = "continuous",
      notes = paste(
        "Screened; not significant. Median 97 mL/min (range 35-304), Table 2.",
        "The categorical renal-function stratification derived from it",
        "(normal / mild / moderate) was also screened and not retained."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator", units = "(binary)", type = "binary",
      notes = "Screened as part of the race covariate; not significant. 2.0% of patients, Table 2."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator", units = "(binary)", type = "binary",
      notes = paste(
        "Screened as part of the race covariate; not significant. 17.0% of",
        "patients (14.3% Japanese, 2.7% other Asian), Table 2. Study 4 was an",
        "all-Japanese dose-escalation study."
      )
    ),
    SMOKING = list(
      description = "Current smoking status indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant. Not tabulated in the main text (ESM Table 2)."
    ),
    CONMED_FULVESTRANT = list(
      description = "Concomitant fulvestrant indicator", units = "(binary)", type = "binary",
      notes = paste(
        "Screened; not significant, unlike the paclitaxel arm which was retained.",
        "16.8% of patients, Table 2."
      )
    ),
    CONMED_CYP3A4_INHIB = list(
      description = "Concomitant CYP3A inhibitor indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A inducer indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant."
    ),
    CONMED_ARA = list(
      description = "Concomitant acid-reducing agent indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant. 31.5% of patients, Sect. 3.1."
    ),
    HEPATIC_IMPAIR = list(
      description = "NCI-ODWG hepatic function category", units = "(category)", type = "categorical",
      notes = paste(
        "Screened; not significant. Normal 67.1%, mild 31.3%, moderate 1.4%,",
        "no severe, Table 2."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 441L,
    n_studies      = 4L,
    n_observations = 3963L,
    age_range      = "27-87 years",
    age_median     = "56 years",
    weight_range   = "32-129 kg",
    weight_median  = "67 kg",
    sex_female_pct = 79.6,
    race_ethnicity = c(
      White = 74.1, Black = 2.0, Asian = 17.0,
      `American Indian or Alaska Native` = 2.9, Other = 2.9, Missing = 0.9
    ),
    disease_state  = "Advanced or metastatic solid tumours (including advanced / metastatic breast cancer)",
    dose_range     = paste(
      "80-800 mg orally twice daily over 21-day and 28-day cycles, as",
      "monotherapy or with paclitaxel or fulvestrant; continuous dosing or one",
      "of two intermittent schedules, 4 days on / 3 days off (4/3, 73.7% of",
      "patients) or 2 days on / 5 days off (2/5, 10.9%)"
    ),
    renal_function = "Normal (CrCL >= 90 mL/min) 58.0%, mild (60-89) 32.9%, moderate (30-59) 8.8%, no severe impairment",
    hepatic_function = "Normal 67.1%, mild 31.3%, moderate 1.4%, no severe impairment",
    co_medication  = "Paclitaxel 20.4% (BEECH only), fulvestrant 16.8% (Study 1 only), acid-reducing agent 31.5%",
    regions        = "Global (Study 1, BEECH and OAK multinational; Study 4 all-Japanese)",
    notes          = paste(
      "Pooled from four phase I / II studies: Study 1 (n = 280), BEECH",
      "(n = 90), Study 4 (n = 41, Japanese) and OAK (n = 30). Baseline",
      "demographics are Fernandez-Teruel 2024 Table 2 (key covariates) and ESM",
      "Table 2 (other covariates). 559 of 4522 samples were excluded, 141 of",
      "them below the 1.00 ng/mL limit of quantification."
    )
  )

  ini({
    # --- Structural disposition -------------------------------------------
    # All parameters are apparent (/F): capivasertib was dosed orally only and
    # no absolute-bioavailability arm is in this dataset, so F is not
    # separately identifiable and is absorbed into CL, V and Q.
    lcl  <- log(62.2);  label("Initial apparent clearance CL0/F, before auto-inhibition (L/h)")                # Fernandez-Teruel 2024 Table 3: CL0/F = 62.2 L/h (RSE 1.58%, bootstrap 95% CI 57.9-66.1)
    lvc  <- log(47.9);  label("Apparent central volume V2/F (L)")                                              # Fernandez-Teruel 2024 Table 3: V2/F = 47.9 L (RSE 1.34%, bootstrap 95% CI 32.4-65.8)
    lvp  <- log(113);   label("Apparent first peripheral volume V3/F (L)")                                     # Fernandez-Teruel 2024 Table 3: V3/F = 113 L (RSE 4.27%, bootstrap 95% CI 79.4-241)
    lq   <- log(2.66);  label("Apparent intercompartmental clearance to peripheral1, Q3/F (L/h)")              # Fernandez-Teruel 2024 Table 3: Q3/F = 2.66 L/h (RSE 1.5%, bootstrap 95% CI 0.941-27)
    lvp2 <- log(94.7);  label("Apparent second peripheral volume V4/F (L)")                                    # Fernandez-Teruel 2024 Table 3: V4/F = 94.7 L (RSE 1.32%, bootstrap 95% CI 52.6-192)
    lq2  <- log(21.8);  label("Apparent intercompartmental clearance to peripheral2, Q4/F (L/h)")              # Fernandez-Teruel 2024 Table 3: Q4/F = 21.8 L/h (RSE 1.31%, bootstrap 95% CI 1-29.1)

    # --- Absorption --------------------------------------------------------
    # Two parallel input routes from the same dose (Fernandez-Teruel 2024
    # Fig. 2): a fraction ffo enters a first-order depot with rate Ka and lag
    # time ALAG1, and the remaining 1 - ffo is delivered as a zero-order input
    # straight into the central compartment over a duration D2.
    lka <- log(0.417);  label("First-order absorption rate constant Ka (1/h)")                                 # Fernandez-Teruel 2024 Table 3: Ka = 0.417 1/h (RSE 1.46%, bootstrap 95% CI 0.278-0.478)
    ld1 <- log(45.1);   label("Duration of the zero-order input into central, D2 (h)")                         # Fernandez-Teruel 2024 Table 3: D2 = 45.1 h (RSE 0.331%, bootstrap 95% CI 42.3-77.7)
    # The paper estimates the first-order fraction on the LOGIT scale
    # (Fernandez-Teruel 2024 Eq. (a)), so `logitffo` -- not `lffo` -- is the
    # correct encoding, and the between-subject random effect is additive on
    # the logit scale as the paper states. logit^-1(1.4) = 0.802, matching the
    # paper's "80% of the dose was absorbed using the first-order mechanism"
    # (Sect. 3.2).
    logitffo <- 1.4;    label("Logit of the fraction of the dose absorbed by the first-order route, F1 (unitless)")  # Fernandez-Teruel 2024 Table 3: Logit F1 = 1.4 (RSE 6.96%, bootstrap 95% CI 1.19-2.75)

    # Two absorption lag times, selected by formulation. Both apply only when
    # the dose was NOT taken after an overnight fast (Fernandez-Teruel 2024
    # Eq. (c) and Sect. 3.2: "There was no Lag when capivasertib was
    # administered orally as a tablet after an overnight fast").
    ltlag_tab <- log(0.212); label("Absorption lag time, tablet, non-overnight-fasted (h)")                    # Fernandez-Teruel 2024 Table 3: Lag1_tab = 0.212 h (RSE 2.11%, bootstrap 95% CI 0.0488-0.314)
    ltlag_cap <- log(0.46);  label("Absorption lag time, capsule, non-overnight-fasted (h)")                   # Fernandez-Teruel 2024 Table 3: Lag1_cap = 0.46 h (RSE 0.292%, bootstrap 95% CI 0.429-0.488)

    # --- Time-dependent auto-inhibition of CL/F ----------------------------
    # Fernandez-Teruel 2024 Eq. (e):
    #   CL/F(t) = CL0/F * (1 - exp(Imax_i) * t^5 / (T50^5 + t^5))
    # The paper's `Imax` is therefore the LOG of the maximal fractional
    # inhibition, not the fraction itself: exp(-1.54) = 0.214, i.e. CL/F falls
    # by at most 21.4% at the 480 mg reference dose. That reading is what
    # reproduces all three of the paper's own reported reductions -- 18%, 22%
    # and 54% at 400, 480 and 800 mg (Sect. 3.5) -- and is the reason the value
    # is stored here on the log scale rather than back-transformed.
    lcl_time_max  <- -1.54;      label("log of the maximal fractional inhibition of CL/F at t >> T50 (unitless); sign applied in model()")  # Fernandez-Teruel 2024 Table 3: Imax = -1.54 (RSE 1.58%, bootstrap 95% CI -1.85 to -1.29)
    lcl_t50       <- log(67.4);  label("T50, time at which half of the maximal inhibition of CL/F is reached (h)")                          # Fernandez-Teruel 2024 Table 3: T50 = 67.4 h (RSE 2.06%, bootstrap 95% CI 54.4-148)
    lcl_time_hill <- fixed(log(5)); label("Hill shape coefficient of the time-on-CL/F function (unitless)")                                 # Fernandez-Teruel 2024 Sect. 3.2: "a Hill parameter fixed at five"; the exponent 5 also appears literally in Eq. (e)

    # --- Covariate effects -------------------------------------------------
    e_wt_cl            <- 0.00585;  label("Linear-deviation coefficient of body weight on CL0/F, per kg above 67 kg (1/kg)")  # Fernandez-Teruel 2024 Table 3: CL0_BBW = 0.00585 (RSE 1.69%, bootstrap 95% CI 0.00289-0.00978); Eq. (d)
    e_wt_logitffo      <- -1.11;    label("Power exponent of (WT / 67) acting on the logit-scale first-order fraction (unitless)")  # Fernandez-Teruel 2024 Table 3: F1_BBW = -1.11 (RSE 20.7%, bootstrap 95% CI -1.56 to -0.618); Eq. (a)
    e_dose_cl_time_max <- -0.00183; label("Linear-deviation coefficient of planned dose on log-Imax, per mg above 480 mg (1/mg)")  # Fernandez-Teruel 2024 Table 3: Imax_dose = -0.00183 (RSE 2.39%, bootstrap 95% CI -0.00241 to -0.00117); Eq. (f)
    e_pacl_cl_time_max <- 1.15;     label("Fractional change in log-Imax with concomitant paclitaxel (unitless)")  # Fernandez-Teruel 2024 Table 3: Imax_pacl = 1.15 (RSE 12.5%, bootstrap 95% CI 1-6.91); Eq. (f)

    # --- Between-subject variability ---------------------------------------
    # Fernandez-Teruel 2024 Table 3 reports every BSV as a percent coefficient
    # of variation, including the one the paper describes as ADDITIVE (Logit
    # F1, Sect. 3.2). A CV can only be read as 100 * sqrt(omega^2) for an
    # additive random effect, so that reporting convention -- rather than the
    # log-normal 100 * sqrt(exp(omega^2) - 1) -- is applied uniformly here:
    # omega = CV / 100 for every eta below.
    etalcl          ~ 0.393^2;    # Fernandez-Teruel 2024 Table 3: CL0/F BSV CV = 39.3% (RSE 5.24%, bootstrap 95% CI 34.5-45)
    etalvc          ~ 1.14^2;     # Fernandez-Teruel 2024 Table 3: V2/F BSV CV = 114% (RSE 7.94%, bootstrap 95% CI 81-133)
    etalcl_time_max ~ 0.706^2;    # Fernandez-Teruel 2024 Table 3: Imax BSV CV = 70.6% (RSE 10.9%, bootstrap 95% CI 53.1-86.4)
    etaltlag        ~ 0.63^2;     # Fernandez-Teruel 2024 Table 3: Lag1 BSV CV = 63% (RSE 10%, bootstrap 95% CI 50.1-74.3)
    etalogitffo     ~ 0.962^2;    # Fernandez-Teruel 2024 Table 3: Logit F1 BSV CV = 96.2% (RSE 16.1%, bootstrap 95% CI 82.5-163); additive on the logit scale per Sect. 3.2
    etalka          ~ fixed(0.15^2);  # Fernandez-Teruel 2024 Table 3 footnote: BSV for Ka and D2 held at 15% CV rather than estimated
    etald1          ~ fixed(0.15^2);  # Fernandez-Teruel 2024 Table 3 footnote: BSV for Ka and D2 held at 15% CV rather than estimated

    # --- Residual unexplained variability ----------------------------------
    # Fitted on log-transformed data with a combined proportional + additive
    # error (Sect. 3.2); rendered here on the linear concentration scale, which
    # is nlmixr2's `prop() + add()` form.
    propSd <- 0.439;  label("Proportional residual error SD (fraction)")        # Fernandez-Teruel 2024 Table 3: RUV proportional CV = 43.9% (RSE 1.39%, bootstrap 95% CI 39.7-48.3)
    addSd  <- 0.504;  label("Additive residual error SD (ng/mL)")               # Fernandez-Teruel 2024 Table 3: RUV additive = 0.504 ug/L (RSE 15.3%, bootstrap 95% CI 0.177-0.721); 1 ug/L = 1 ng/mL
  })

  model({
    # --- 1. Absorption -----------------------------------------------------
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    # Fraction absorbed by the first-order route, on the logit scale, with a
    # power effect of body weight and an additive between-subject random
    # effect. Fernandez-Teruel 2024 Eq. (a):
    #   F1i = exp(LogitF1 * (BBW/67)^F1_BBW + eta) /
    #         (1 + exp(LogitF1 * (BBW/67)^F1_BBW + eta))
    # The eta is kept on its own line so rxode2 sees a mu-referenced term.
    logitffo_ind <- logitffo * (WT / 67)^e_wt_logitffo + etalogitffo
    ffo <- 1 / (1 + exp(-logitffo_ind))

    # Absorption lag time. Formulation selects the magnitude and the overnight
    # fast switches it off entirely. Fernandez-Teruel 2024 Eq. (c):
    #   ALAG1i = (Lag1_tab*(1-CAP) + Lag1_cap*CAP) * (1 - FASTED) * exp(eta)
    tlag <- (exp(ltlag_tab) * (1 - FORM_CAPSULE) + exp(ltlag_cap) * FORM_CAPSULE) *
      (1 - FASTED_STRICT) * exp(etaltlag)

    # --- 2. Disposition ----------------------------------------------------
    # Fernandez-Teruel 2024 Eq. (d): CL0/Fi = CL0/F * (1 + (BBW - 67)*CL0_BBW) * exp(eta)
    cl_base <- exp(lcl + etalcl) * (1 + (WT - 67) * e_wt_cl)
    vc      <- exp(lvc + etalvc)
    vp      <- exp(lvp)
    q       <- exp(lq)
    vp2     <- exp(lvp2)
    q2      <- exp(lq2)

    # --- 3. Time-dependent auto-inhibition of CL/F -------------------------
    # Fernandez-Teruel 2024 Eq. (f): the dose and paclitaxel effects are
    # MULTIPLICATIVE on the log-scale magnitude, and so is the random effect:
    #   Imax_i = Imax * (1 + (DOSE - 480)*Imax_dose) * (1 + PACL*Imax_pacl) * exp(eta)
    # Since Imax < 0, a multiplier > 1 deepens the inhibition and a multiplier
    # < 1 shallows it.
    lcl_time_max_ind <- lcl_time_max *
      (1 + (DOSE_CAPIVASERTIB_MG - 480) * e_dose_cl_time_max) *
      (1 + CONMED_PACLITAXEL * e_pacl_cl_time_max) *
      exp(etalcl_time_max)
    cl_time_max_i <- -exp(lcl_time_max_ind)
    cl_t50        <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)

    # Fernandez-Teruel 2024 Eq. (e). `t` is time since the first dose.
    cl <- cl_base *
      (1 + cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # --- 4. Micro-constants and ODEs ---------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # --- 5. Parallel input routing -----------------------------------------
    # The same dose is recorded twice: a fraction `ffo` enters `depot` as a
    # lagged first-order input, and the complementary fraction 1 - ffo
    # (Fernandez-Teruel 2024 Eq. (b), F2i = 1 - F1i) enters `central` as a
    # zero-order input of duration d1. The central dose record must request a
    # modelled duration (rate = -2 in the event table).
    f(depot)     <- ffo
    alag(depot)  <- tlag
    f(central)   <- 1 - ffo
    dur(central) <- d1

    # --- 6. Observation ----------------------------------------------------
    # Dose is in mg and vc in L, so central/vc is mg/L; x 1000 gives ng/mL
    # (= ug/L), the unit Fernandez-Teruel 2024 reports concentrations in.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
