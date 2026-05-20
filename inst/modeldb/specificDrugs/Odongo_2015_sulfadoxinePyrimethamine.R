Odongo_2015_sulfadoxinePyrimethamine <- function() {
  description <- paste(
    "Joint popPK model for the antimalarial fixed-dose combination of",
    "sulfadoxine (1500 mg) and pyrimethamine (75 mg) administered as a",
    "single oral dose for intermittent preventive treatment of malaria",
    "during pregnancy (IPTp) in 34 non-pregnant and 87 pregnant Ugandan",
    "women dosed in the second trimester, of whom 78 were redosed in",
    "the third trimester (Odongo 2015). Each drug is described by a",
    "two-compartment model with first-order absorption and an",
    "absorption lag time, with bioavailability fixed at 1. Covariates",
    "on apparent CL/F (additive in L/h): pregnancy status (both",
    "drugs), serum albumin (sulfadoxine only), and subject age",
    "(pyrimethamine only). Covariates on apparent central volume V2/F",
    "(exponential per-unit): gestational age at dose (both drugs) and",
    "body weight (pyrimethamine only). Inter-individual variability is",
    "log-normal and is not estimated on V2/F or V3/F for sulfadoxine,",
    "nor on Q/F for pyrimethamine, in line with the paper's",
    "over-parameterisation control."
  )
  reference <- paste(
    "Odongo CO, Bisaso KR, Ntale M, Odia G, Ojara FW, Byamugisha J,",
    "Mukonzo JK, Obua C.",
    "Trimester-Specific Population Pharmacokinetics and Other",
    "Correlates of Variability in Sulphadoxine-Pyrimethamine",
    "Disposition Among Ugandan Pregnant Women.",
    "Drugs R D. 2015 Dec;15(4):351-362.",
    "doi:10.1007/s40268-015-0110-z.",
    sep = " "
  )
  vignette <- "Odongo_2015_sulfadoxinePyrimethamine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L (= ug/mL) for sulfadoxine plasma; ng/mL for pyrimethamine plasma"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at dosing",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per dosing occasion (re-collected when women were",
        "redosed in the third trimester). Pyrimethamine apparent V2/F",
        "is scaled exponentially as V2 = V2_TV * exp(0.0084 * (WT -",
        "60)) with reference WT = 60 kg (the population median, also",
        "the second-trimester median; Table 1 of Odongo 2015). Body",
        "weight was tested as a covariate on sulfadoxine V2/F but did",
        "not reach statistical significance and is not retained for",
        "sulfadoxine."
      ),
      source_name        = "Body weight (kg)"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per dosing occasion. Pyrimethamine apparent CL/F",
        "is shifted additively by 0.016 L/h per year of age, centred",
        "at the population-median reference 23 years (cohort mean",
        "across 199 dosing occasions; Table 1 of Odongo 2015 reports",
        "23.7 / 22.8 / 23.5 years across non-pregnant / T2 / T3).",
        "Age was tested as a covariate on sulfadoxine CL/F but was not",
        "retained in the final model."
      ),
      source_name        = "Age (years)"
    ),
    ALB = list(
      description        = "Serum albumin concentration at dosing",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per dosing occasion. Sulfadoxine apparent CL/F is",
        "shifted additively by 0.013 L/h per unit decrease in serum",
        "albumin (Odongo 2015 Discussion: 'a unit decrease in the",
        "serum albumin level was associated with an increase in CL/F",
        "of 0.013 L/h'), centred at the non-pregnant reference 44.6",
        "g/L (Table 1 mean for the non-pregnant cohort). Albumin units",
        "are recorded as g/L based on the biological plausibility of",
        "the reported Table 1 values (44.6 / 37.4 / 34.1) which only",
        "fall in the physiologic range when interpreted as g/L; the",
        "Table 1 header 'g/dL' is treated as a publication typo (see",
        "vignette Errata). The covariate had less than 50% inclusion",
        "frequency at the bootstrap stage but was retained by the",
        "authors for biologic plausibility; the typical-value",
        "predictions produced by this term during pregnancy disagree",
        "with the per-trimester Bayesian post-hoc estimates in Table",
        "3 of the paper, also flagged in vignette Errata."
      ),
      source_name        = "Serum albumin (g/L per Table 1 plausibility; reported as g/dL)"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester); 0 = non-pregnant.",
        "Both sulfadoxine and pyrimethamine apparent CL/F are shifted",
        "additively (in L/h) by a pregnancy step: +0.0284 L/h for",
        "sulfadoxine and +0.319 L/h for pyrimethamine (Table 2,",
        "Pregnancy-CL rows). Reference category is the non-pregnant",
        "state."
      ),
      source_name        = "Pregnancy status"
    ),
    GA = list(
      description        = "Gestational age at the time of dosing",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject across the second- and",
        "third-trimester visits (per-occasion gestational age at",
        "dose, NOT gestational age at birth). For non-pregnant",
        "subjects GA is set to 0 in the dataset so the exponential",
        "centring (GA - 20) cleanly produces the non-pregnant typical",
        "V2/F. Both sulfadoxine V2/F (coefficient 0.0093 / week) and",
        "pyrimethamine V2/F (coefficient 0.0079 / week) are scaled",
        "exponentially with reference GA = 20 weeks (the",
        "second-trimester median; Table 1 of Odongo 2015 reports",
        "median GA at dosing of 20 weeks in T2 and 28 weeks in T3).",
        "Semantic note: the nlmixr2lib canonical 'GA' is documented",
        "for neonatal / paediatric models as gestational age at",
        "birth (time-fixed); the same column name is used here for a",
        "maternal-pregnancy popPK model where GA is per-occasion. The",
        "unit (weeks) and biological concept (pregnancy duration",
        "from menstrual start) are identical; only the time of",
        "recording differs."
      ),
      source_name        = "Gestational age (weeks)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 199,
    n_subjects_unique = 121,
    n_studies      = 1,
    age_range      = paste(
      "Adult women of reproductive age (means 23.7 / 22.8 / 23.5",
      "years across non-pregnant / T2 / T3 cohorts; Table 1 of",
      "Odongo 2015)."
    ),
    weight_range   = paste(
      "Adult women (per-group median 58.0 / 60.0 / 63.5 kg across",
      "non-pregnant / T2 / T3; cohort median 60 kg, the second-",
      "trimester median; Table 1 of Odongo 2015)."
    ),
    sex_female_pct = 100,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy Ugandan women dosed as part of an intermittent",
      "preventive treatment in pregnancy (IPTp) study. Inclusion:",
      "pregnancy of at least 16 weeks (median GA 20 weeks in T2 and",
      "28 weeks in T3) or non-pregnant; HIV negative; no SP use in",
      "the current pregnancy; no SP allergy; no chronic illness",
      "directly or indirectly related to pregnancy. Median",
      "gestational age at dosing was 20 weeks in T2 and 28 weeks in",
      "T3; the same women were sequentially redosed in T3 after a",
      "wash-out of 2-6 weeks following the T2 follow-up window."
    ),
    dose_range     = paste(
      "Single oral dose of three fixed-dose combination tablets",
      "(Malaren brand; Rene Industries, Kampala, Uganda) containing",
      "a total of 1500 mg sulfadoxine and 75 mg pyrimethamine,",
      "swallowed whole with plain water under midwife supervision",
      "after an overnight fast. Subjects in the pregnant cohort",
      "received a second matched dose during the third trimester."
    ),
    regions        = "East Africa (Uganda; Mulago Hospital, Kampala)",
    notes          = paste(
      "Demographics and dataset composition from Odongo 2015 Table",
      "1 and Section 2.5.1. 34 non-pregnant women contributed 172",
      "sulfadoxine and 172 pyrimethamine concentration records",
      "above time zero; 87 pregnant women dosed in T2 contributed",
      "425 / 418 records respectively; 78 of those women were",
      "redosed in T3, contributing 378 / 372 records and treated as",
      "distinct individuals in the analysis (199 dosing occasions",
      "total). Sampling times ranged from 0.5 h to 42 days postdose",
      "(median ~5 observations per dosing occasion). Concentrations",
      "were quantified in plasma by HPLC-UV (LLOQ 25 umol/L for",
      "sulfadoxine and 40 ng/mL for pyrimethamine). The fit was done",
      "in NONMEM v7.2 with FOCE-I using log-transformed",
      "concentrations and an additive residual error on the log",
      "scale (log-normal in linear space)."
    )
  )

  ini({
    # ============================================================
    # Sulfadoxine structural parameters (parent drug)
    # ----- Odongo 2015 Table 2, "Sulphadoxine Final (RSE %)" column
    # All apparent volumes and flow rates are in plasma units, with
    # bioavailability F1 = 1 by structural assumption.
    # ============================================================
    lka     <- log(0.664)
    label("Sulfadoxine first-order absorption rate, ka (1/h)")            # Table 2: ka_TV = 0.664 /h (RSE 18.5%)
    lcl     <- log(0.0059)
    label("Sulfadoxine apparent CL/F, typical non-pregnant baseline at ALB = 44.6 g/L (L/h)") # Table 2: CL/F_TV = 0.0059 L/h (RSE 27.8%)
    lvc     <- log(10.74)
    label("Sulfadoxine apparent central volume V2/F at reference GA = 20 weeks (L)") # Table 2: V2/F_TV = 10.74 L (RSE 2.8%)
    lq      <- log(0.029)
    label("Sulfadoxine apparent inter-compartmental CL (L/h)")            # Table 2: Q_TV = 0.029 L/h (RSE 17.4%)
    lvp     <- log(161.71)
    label("Sulfadoxine apparent peripheral volume V3/F (L)")              # Table 2: V3/F_TV = 161.71 L (RSE 33.8%)
    lalag   <- log(0.371)
    label("Sulfadoxine absorption lag time, ALAG (h)")                    # Table 2: ALAG_TV = 0.371 h (RSE 11.8%)
    lfdepot <- fixed(log(1))
    label("Sulfadoxine relative bioavailability F (unitless, FIXED at 1)") # Section 3.1: "bioavailability was assumed to be equal to 1"

    # ============================================================
    # Pyrimethamine structural parameters
    # ----- Odongo 2015 Table 2, "Pyrimethamine Final (RSE %)" column
    # ============================================================
    lka_pyra     <- log(1.216)
    label("Pyrimethamine first-order absorption rate, ka (1/h)")              # Table 2: ka_TV = 1.216 /h (RSE 5.7%)
    lcl_pyra     <- log(0.545)
    label("Pyrimethamine apparent CL/F, typical non-pregnant baseline at AGE = 23 y (L/h)") # Table 2: CL/F_TV = 0.545 L/h (RSE 5.0%)
    lvc_pyra     <- log(153.915)
    label("Pyrimethamine apparent central volume V2/F at reference GA = 20 weeks, WT = 60 kg (L)") # Table 2: V2/F_TV = 153.915 L (RSE 2.2%)
    lq_pyra      <- log(0.297)
    label("Pyrimethamine apparent inter-compartmental CL (L/h)")              # Table 2: Q_TV = 0.297 L/h (RSE 30.9%)
    lvp_pyra     <- log(51.224)
    label("Pyrimethamine apparent peripheral volume V3/F (L)")                # Table 2: V3/F_TV = 51.224 L (RSE 16.1%)
    lalag_pyra   <- log(0.394)
    label("Pyrimethamine absorption lag time, ALAG (h)")                      # Table 2: ALAG_TV = 0.394 h (RSE 1.3%)
    lfdepot_pyra <- fixed(log(1))
    label("Pyrimethamine relative bioavailability F (unitless, FIXED at 1)")  # Section 3.2: same structural assumption as sulfadoxine

    # ============================================================
    # Covariate effects on apparent CL/F (additive in L/h)
    # Sign convention: a positive coefficient increases CL/F when the
    # covariate moves toward the direction stated by the Odongo 2015
    # Table 2 footnote (pregnancy ON, age UP, albumin DOWN).
    # ============================================================
    e_preg_cl <- 0.0284
    label("Sulfadoxine additive CL/F shift for pregnancy (L/h)")              # Table 2: Pregnancy-CL_b = 0.0284 (RSE 16.3%); footnote b "Increase in CL with pregnancy"
    e_alb_cl  <- 0.013
    label("Sulfadoxine additive CL/F shift per g/L decrease in serum albumin (L/h)") # Table 2: Albumin-CL_a = 0.013 (RSE 11.5%); footnote a "Increase in CL with unit decrease in albumin level"

    e_preg_cl_pyra <- 0.319
    label("Pyrimethamine additive CL/F shift for pregnancy (L/h)")            # Table 2: Pregnancy-CL_b = 0.319 (RSE 15.8%); footnote b
    e_age_cl_pyra  <- 0.016
    label("Pyrimethamine additive CL/F shift per year of age (L/h)")          # Table 2: Age-CL = 0.016 (RSE 36.2%)

    # ============================================================
    # Covariate effects on apparent central V2/F (exponential per
    # unit covariate, centred at the reference). Form is
    # V2 = V2_TV * exp(coef * (cov - cov_ref)), per Table 2 footnotes
    # c and d ("Exponential increase in V2 with [...] increase in [GA / WT]").
    # ============================================================
    e_ga_vc <- 0.0093
    label("Sulfadoxine exponential coefficient on V2/F per week of gestational age (1/week)") # Table 2: Gestation-V2_c = 0.0093 (RSE 14.1%)

    e_ga_vc_pyra <- 0.0079
    label("Pyrimethamine exponential coefficient on V2/F per week of gestational age (1/week)") # Table 2: Gestation-V2_c = 0.0079 (RSE 23.1%)
    e_wt_vc_pyra <- 0.0084
    label("Pyrimethamine exponential coefficient on V2/F per kg of body weight (1/kg)")        # Table 2: Body weight-V2_d = 0.0084 (RSE 23.7%)

    # ============================================================
    # Inter-individual variability. omega^2 = log(1 + CV^2) for the
    # log-normal eta on each log-transformed parameter; CV% values
    # come from Table 2 ("Final (RSE %)" column, IIV rows). The
    # paper held V2/F, V3/F, and ALAG IIV for sulfadoxine fixed at
    # zero to avoid over-parameterisation; pyrimethamine omits Q/F
    # IIV for the same reason. ALAG IIV is not estimated for either
    # drug (Section 3.1 / 3.2). Those etas are therefore not
    # declared in ini().
    # ============================================================
    etalka      ~ log(1 + 1.026^2)
    # Table 2 sulfa: IIV_KA = 102.6% CV (RSE 11.3%) -> omega^2 = log(1 + 1.026^2)
    etalcl      ~ log(1 + 0.446^2)
    # Table 2 sulfa: IIV_CL = 44.6% CV (RSE 17.7%) -> omega^2 = log(1 + 0.446^2)
    etalq       ~ log(1 + 0.523^2)
    # Table 2 sulfa: IIV_Q = 52.3% CV (RSE 17.5%) -> omega^2 = log(1 + 0.523^2)

    etalka_pyra ~ log(1 + 1.209^2)
    # Table 2 pyra: IIV_KA = 120.9% CV (RSE 10.3%) -> omega^2 = log(1 + 1.209^2)
    etalcl_pyra ~ log(1 + 0.305^2)
    # Table 2 pyra: IIV_CL = 30.5% CV (RSE 13.3%) -> omega^2 = log(1 + 0.305^2)
    etalvc_pyra ~ log(1 + 0.081^2)
    # Table 2 pyra: IIV_V2 = 8.1% CV (RSE 62.9%) -> omega^2 = log(1 + 0.081^2)
    etalvp_pyra ~ log(1 + 1.093^2)
    # Table 2 pyra: IIV_V3 = 109.3% CV (RSE 12.9%) -> omega^2 = log(1 + 1.093^2)

    # ============================================================
    # Residual error - the paper fitted log-transformed concentration
    # data with an additive residual on the log scale, which is the
    # canonical log-normal residual on the linear-scale prediction
    # (Section 2.5.1). nlmixr2's lnorm() captures exactly this
    # structure; the reported "Residual (CV %)" rows of Table 2 are
    # the SD of the log-additive residual on log scale, equivalent to
    # the approximate linear-space CV.
    # ============================================================
    expSd      <- 0.331
    label("Sulfadoxine log-normal residual SD (unitless on log scale)")       # Table 2 sulfa: Residual = 33.1% CV (RSE 8.5%)
    expSd_pyra <- 0.272
    label("Pyrimethamine log-normal residual SD (unitless on log scale)")     # Table 2 pyra: Residual = 27.2% CV (RSE 7.5%)
  })

  model({
    # ------------------------------------------------------------
    # Covariate reference values. All chosen to match the typical
    # individual whose covariates produce the Table 2 TV values
    # (back-calculated against the per-trimester Bayesian post-hoc
    # estimates in Table 3 of the paper).
    # ------------------------------------------------------------
    ga_ref  <- 20    # weeks; second-trimester median, Table 1
    wt_ref  <- 60    # kg; cohort median (also T2 median), Table 1
    age_ref <- 23    # years; cohort-weighted mean age, Table 1
    alb_ref <- 44.6  # g/L; non-pregnant mean, Table 1

    # ------------------------------------------------------------
    # Sulfadoxine individual PK parameters.
    # Additive covariates on CL/F are applied to the typical value
    # exp(lcl) before multiplying by the log-normal eta on CL.
    # GA enters V2/F exponentially; non-pregnant subjects have GA = 0
    # in the dataset so the exponential evaluates to exp(0.0093 *
    # (0 - 20)) ~ 0.83, recovering V2_nonpreg = 8.92 L (Table 3).
    # ------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl_tv <- exp(lcl) +
             e_preg_cl * PREG +
             e_alb_cl  * (alb_ref - ALB)
    cl <- cl_tv * exp(etalcl)
    vc <- exp(lvc) * exp(e_ga_vc * (GA - ga_ref))
    q  <- exp(lq + etalq)
    vp <- exp(lvp)
    alag <- exp(lalag)

    # ------------------------------------------------------------
    # Pyrimethamine individual PK parameters.
    # ------------------------------------------------------------
    ka_pyra <- exp(lka_pyra + etalka_pyra)
    cl_tv_pyra <- exp(lcl_pyra) +
                  e_preg_cl_pyra * PREG +
                  e_age_cl_pyra  * (AGE - age_ref)
    cl_pyra <- cl_tv_pyra * exp(etalcl_pyra)
    vc_pyra <- exp(lvc_pyra + etalvc_pyra) *
               exp(e_ga_vc_pyra * (GA - ga_ref)) *
               exp(e_wt_vc_pyra * (WT - wt_ref))
    q_pyra  <- exp(lq_pyra)
    vp_pyra <- exp(lvp_pyra + etalvp_pyra)
    alag_pyra <- exp(lalag_pyra)

    # ------------------------------------------------------------
    # Sulfadoxine 2-compartment disposition with first-order
    # absorption from a `depot` compartment and an absorption-lag
    # time. Compartment names use the canonical parent-drug naming
    # (depot, central, peripheral1, Cc).
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl + q) * central / vc +
                          q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp

    f(depot)    <- exp(lfdepot)
    alag(depot) <- alag

    # ------------------------------------------------------------
    # Pyrimethamine 2-compartment disposition with first-order
    # absorption from a separate `depot_pyra` compartment and its
    # own absorption-lag time. Compartment names follow the
    # parent + metabolite-suffix convention with the `_pyra` token.
    # ------------------------------------------------------------
    d/dt(depot_pyra)       <- -ka_pyra * depot_pyra
    d/dt(central_pyra)     <-  ka_pyra * depot_pyra -
                               (cl_pyra + q_pyra) * central_pyra / vc_pyra +
                               q_pyra * peripheral1_pyra / vp_pyra
    d/dt(peripheral1_pyra) <-  q_pyra * central_pyra / vc_pyra -
                               q_pyra * peripheral1_pyra / vp_pyra

    f(depot_pyra)    <- exp(lfdepot_pyra)
    alag(depot_pyra) <- alag_pyra

    # ------------------------------------------------------------
    # Plasma concentration outputs. Sulfadoxine in mg/L (= ug/mL);
    # pyrimethamine in ng/mL (mg/L * 1000), matching the assay-
    # reporting unit (LLOQ 40 ng/mL) and Table 3's pyrimethamine
    # AUC reporting unit.
    # ------------------------------------------------------------
    Cc      <- central      / vc
    Cc_pyra <- (central_pyra / vc_pyra) * 1000

    Cc      ~ lnorm(expSd)
    Cc_pyra ~ lnorm(expSd_pyra)
  })
}
