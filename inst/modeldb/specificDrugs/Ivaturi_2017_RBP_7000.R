# Joint population PK / PANSS PD model for the once-monthly long-acting
# subcutaneous Atrigel formulation of risperidone (RBP-7000) in 354 adults
# with acute schizophrenia (Ivaturi 2017, Br J Clin Pharmacol 83(7):1476-1498;
# doi:10.1111/bcp.13246).

Ivaturi_2017_RBP_7000 <- function() {
  description <- paste(
    "Integrated population pharmacokinetic / PANSS pharmacodynamic model for",
    "the once-monthly long-acting subcutaneous Atrigel formulation of",
    "risperidone (RBP-7000, Indivior) in 337 adults with acute schizophrenia",
    "treated with two SC injections (90 mg or 120 mg) 28 days apart in a Phase",
    "3 registration trial (NCT02109562). The PK sub-model is the empirical",
    "dual-absorption structure inherited from the upstream RBP-7000 SAD and",
    "MAD studies (Gomeni 2013; Laffont 2014, 2015): a fast first-order",
    "absorption rate ka1 from the SC depot to the risperidone central",
    "compartment captures the rapid release from the injection site, while a",
    "5-compartment transit chain with rate constant ktr feeds a slow",
    "first-order absorption rate ka2 from the terminal transit compartment",
    "into central, mimicking the slow sustained release from the solidified",
    "ATRIGEL implant. Systemically available risperidone is distributed to a",
    "single peripheral compartment (rate constants krrp and krpr), eliminated",
    "by non-metabolite routes (krel), and partly converted to its equipotent",
    "9-hydroxyrisperidone metabolite at rate kr9; the metabolite is described",
    "by a one-compartment model with first-order elimination (k9el) and an",
    "apparent volume of distribution constrained equal to the central volume",
    "of the parent V because VM was not identifiable. CYP2D6 intermediate and",
    "poor metabolizers (vs the extensive / inconclusive reference) have 76",
    "and 94 percent lower metabolite formation rate, respectively. Plasma",
    "concentrations of risperidone and 9-OH-risperidone are converted to",
    "total active moiety in risperidone-equivalent units by AM = [risperidone]",
    "+ [9-OH-risperidone] * 410/426 (molecular-weight correction). The PANSS",
    "PD sub-model combines a Weibull-shaped placebo response (PMAX, TPROG,",
    "POW), an additive linear-drift term (DRIFT) that captures the",
    "improvement-then-worsening pattern observed in some individuals, and an",
    "Emax model relating total active moiety to relative PANSS decrease, with",
    "drug and placebo effects entering additively per Predicted PANSS = BSL *",
    "(1 - PMAX * (1 - exp(-(T/TPROG)^POW)) - Emax * AM / (EC50 + AM)) + DRIFT *",
    "T (T in weeks). The proportional-odds CGI-S sub-model of Table 4 is",
    "documented in the vignette but not implemented here because rxode2's",
    "additive residual / d/dt() ODE pipeline does not natively express",
    "ordinal-logistic observation likelihoods; see the vignette Assumptions",
    "and deviations.",
    sep = " "
  )
  reference <- paste(
    "Ivaturi V, Gopalakrishnan M, Gobburu JVS, Zhang W, Liu Y, Heidbreder C,",
    "Laffont CM (2017). Exposure-response analysis after subcutaneous",
    "administration of RBP-7000, a once-a-month long-acting Atrigel",
    "formulation of risperidone. British Journal of Clinical Pharmacology",
    "83(7):1476-1498. doi:10.1111/bcp.13246.",
    sep = " "
  )
  vignette <- "Ivaturi_2017_RBP_7000"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, intermediate, or inconclusive metabolizer; both CYP2D6_PM and CYP2D6_IM = 0 indicates the pooled EM + Inconclusive reference)",
      notes              = paste(
        "1 = subject is a CYP2D6 poor metabolizer (genotype encoding no",
        "functional enzyme activity), 0 otherwise. Paired with CYP2D6_IM to",
        "encode the three-level PM / IM / EM phenotype following the Ivaturi",
        "2017 Methods reference-category choice: 'Extensive Metabolizers (EM)",
        "and Inconclusive = reference, Poor Metabolizers (PM), Intermediate",
        "Metabolizers (IM)'. Both indicators 0 maps to the pooled EM +",
        "Inconclusive reference. Multiplicative effect on the metabolite",
        "formation rate kr9: kr9_PM = kr9_typical * (1 + e_cyp2d6_pm_kmet)",
        "with e_cyp2d6_pm_kmet = -0.94 (Table 2; 94% lower metabolite",
        "formation than EM / Inconclusive). Phenotype distribution in this",
        "Phase 3 cohort: 82-88% EM, 3.6-7.1% IM, 0.9-2.6% PM, 5.2-7.1%",
        "Inconclusive (Table 1)."
      ),
      source_name        = "CYP2D6 phenotype (Table 1; from CYP2D6 genotyping in pharmacogenetic sub-study)"
    ),
    CYP2D6_IM = list(
      description        = "CYP2D6 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, poor, or inconclusive metabolizer; both CYP2D6_PM and CYP2D6_IM = 0 indicates the pooled EM + Inconclusive reference)",
      notes              = paste(
        "1 = subject is a CYP2D6 intermediate metabolizer (genotype encoding",
        "reduced enzyme activity), 0 otherwise. Paired with CYP2D6_PM with",
        "the pooled EM + Inconclusive reference (both indicators 0).",
        "Multiplicative effect on the metabolite formation rate kr9: kr9_IM =",
        "kr9_typical * (1 + e_cyp2d6_im_kmet) with e_cyp2d6_im_kmet = -0.76",
        "(Table 2; 76% lower metabolite formation than EM / Inconclusive)."
      ),
      source_name        = "CYP2D6 phenotype (Table 1; from CYP2D6 genotyping in pharmacogenetic sub-study)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the full PK and PK/PD covariate analyses but not retained in the final model (Results, 'No covariates were identified' on the absorption / disposition / PD parameters)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained (Results)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on absorption rate constants (significant in the upstream Gomeni 2013 and Laffont 2014 papers; lower ka1 in subjects with higher BMI was attributed to fat-abdominal-tissue effects on absorption of the lipophilic risperidone) but not retained in the Ivaturi 2017 Phase 3 final model (Results)."
    ),
    WTH = list(
      description = "Waist-to-hip ratio",
      units       = "unitless",
      type        = "continuous",
      notes       = "Screened on absorption rate constants but not retained (Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a liver-function predictor but not retained (Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a liver-function predictor but not retained (Results)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as a renal-function predictor but not retained (Results)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (male = reference per Methods) but not retained on either PK or PK/PD parameters (Results)."
    ),
    RACE_BLACK = list(
      description = "Black or African-American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was tested in two encodings: graphical EBE-vs-covariate plots",
        "showed a correlation between race and the absorption rate constants",
        "ka1 and ka2 (so race was carried forward into NONMEM covariate",
        "testing) but the formal stepwise selection rejected race on the",
        "absorption rates (Results, 'Race was also not significant'). The",
        "paper's reference category was Black / African-American (the",
        "majority of the cohort, 70-75%; Table 1)."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 337L,
    n_studies       = 1L,
    age_range       = "approximately 18-65 years; mean 41-43 years across treatment arms (Table 1)",
    age_median      = "approximately 41 years (means 42.76, 40.45, 40.41 years for placebo, 90 mg, 120 mg arms)",
    weight_range    = "mean 88.5-92.6 kg; SD 18.9-22.9 kg across arms (Table 1)",
    weight_median   = "approximately 90 kg (mean values reported, not medians)",
    sex_female_pct  = 23.4,
    race_ethnicity  = c(
      `Black or African-American` = 72.1,
      Other                       = 27.9
    ),
    ethnicity       = "91.5% non-Hispanic / non-Latino, 8.3% Hispanic / Latino (Table 1)",
    disease_state   = paste(
      "Adults with acute schizophrenia, acute psychotic state, or relapse",
      "with acute schizophrenic symptoms; baseline PANSS 80-120 and at least",
      "two of four positive subscale items scoring >4 at screening (Methods).",
      "Mean baseline PANSS 94 across arms (Table 3)."
    ),
    dose_range      = paste(
      "Two subcutaneous injections of RBP-7000 (90 mg or 120 mg)",
      "administered 28 days apart on Day 1 and Day 29 (Methods)."
    ),
    cyp2d6_distribution = paste(
      "Table 1: 82.1-88.2% extensive, 3.6-7.1% intermediate, 0.9-2.6% poor,",
      "5.2-7.1% inconclusive, ~1% missing across treatment arms."
    ),
    regions         = "USA (35 clinical sites)",
    notes           = paste(
      "ITT population n = 337: 112 placebo, 111 90 mg, 114 120 mg. 234",
      "subjects in the RBP-7000 treatment arms contributed 3154 PK samples",
      "(1577 risperidone + 1577 9-OH-risperidone) over the 8-week study;",
      "1571 PANSS and 1549 CGI-S assessments were also collected (Day 1,",
      "Day 15, Day 29, Day 43, Day 57). NONMEM 7.3 with FOCE-I for the PK",
      "and continuous PANSS models, LAPLACE for the proportional-odds CGI-S",
      "model. Concomitant medications and clinical study site were not",
      "tested as covariates due to high missingness / low per-site counts",
      "(Results)."
    )
  )

  ini({
    # ------------------- Pharmacokinetic structural parameters --------------
    # Rate constants on the log scale. Source values come from Table A2 of
    # the Appendix (Phase-3 column) which reports the same final-model
    # estimates as Table 2 but to three significant figures; Table 2 rounds
    # to two figures. Time unit is the hour throughout this section.
    lka1     <- log(0.005); label("Log fast first-order absorption rate ka1 (1/h)")                       # Table 2 / A2: ka1 = 0.005 1/h (95% CI 0.004-0.007; rapid release from SC injection site)
    lka2     <- log(0.016); label("Log slow absorption rate ka2 from terminal transit to central (1/h)")  # Table 2 / A2: ka2 = 0.016 1/h (95% CI 0.010-0.022; ATRIGEL sustained release path)
    lktr     <- log(0.023); label("Log transit rate constant ktr (1/h)")                                  # Table 2 / A2: ktr = 0.023 1/h (95% CI 0.021-0.025; 5-compartment transit chain)
    lkel     <- log(0.043); label("Log risperidone elimination rate krel via non-metabolite routes (1/h)") # Table A2 (Phase-3): krel = 0.043 1/h; Table 2 rounded value 0.04 1/h
    lkmet    <- log(0.221); label("Log metabolite formation rate kr9 from risperidone to 9-OH-risperidone (1/h)") # Table A2 (Phase-3): kr9 = 0.221 1/h
    lkel_9oh <- log(0.069); label("Log 9-OH-risperidone elimination rate k9el (1/h)")                     # Table A2 (Phase-3): k9el = 0.069 1/h
    lk12     <- log(0.841); label("Log central-to-peripheral rate constant krrp (1/h)")                   # Table A2 (Phase-3): krrp = 0.841 1/h
    lk21     <- log(0.006); label("Log peripheral-to-central rate constant krpr (1/h)")                   # Table A2 (Phase-3): krpr = 0.006 1/h
    lvc      <- log(129);   label("Log apparent central volume of distribution V (L); shared between risperidone and 9-OH-risperidone (VM identifiability constraint)") # Table A2 (Phase-3): V = 129 L

    # ------------------- CYP2D6 covariate effects ---------------------------
    # Multiplicative proportional shifts on the metabolite formation rate
    # kr9 relative to the pooled extensive-metabolizer + inconclusive
    # reference (Methods 'Covariate model' for categorical covariates:
    # P_TV = P_TV * (1 + prop_factor) if Cov = 1).
    e_cyp2d6_im_kmet <- -0.76; label("Proportional CYP2D6 intermediate-metabolizer shift on kr9 (unitless)") # Table 2: CYP2D6 Intermediate effect on metabolite formation = -0.76 (95% CI -0.85, -0.56)
    e_cyp2d6_pm_kmet <- -0.94; label("Proportional CYP2D6 poor-metabolizer shift on kr9 (unitless)")          # Table 2: CYP2D6 Poor effect on metabolite formation = -0.94 (95% CI -0.96, -0.88)

    # ------------------- Inter-individual variability on PK -----------------
    # IIV on the log-transformed parameter; omega^2 = log(CV^2 + 1) using
    # the NONMEM convention because the paper reports IIV on each PK
    # parameter as the percent-CV calculated as sqrt(exp(omega^2) - 1)
    # (Methods, Equation for exponential models). No IIV is reported on
    # krel.
    etalka1     ~ 0.161  # Table 2: BSV(ka1)   = 41.8% CV; omega^2 = log(1 + 0.418^2) = 0.161
    etalka2     ~ 0.098  # Table 2: BSV(ka2)   = 32.1% CV; omega^2 = log(1 + 0.321^2) = 0.098
    etalktr     ~ 0.165  # Table 2: BSV(ktr)   = 42.3% CV; omega^2 = log(1 + 0.423^2) = 0.165
    etalkmet    ~ 0.218  # Table 2: BSV(kr9)   = 49.3% CV; omega^2 = log(1 + 0.493^2) = 0.218
    etalkel_9oh ~ 0.035  # Table 2: BSV(k9el)  = 18.8% CV; omega^2 = log(1 + 0.188^2) = 0.035
    etalk12     ~ 0.185  # Table 2: BSV(krrp)  = 45.1% CV; omega^2 = log(1 + 0.451^2) = 0.185
    etalk21     ~ 0.379  # Table 2: BSV(krpr)  = 67.9% CV; omega^2 = log(1 + 0.679^2) = 0.379
    etalvc      ~ 0.139  # Table 2: BSV(V)     = 38.6% CV; omega^2 = log(1 + 0.386^2) = 0.139

    # ------------------- Residual error on plasma concentrations ------------
    # Combined additive + proportional error model common to risperidone
    # and 9-OH-risperidone (Methods, Equation 'PK Model' and Table 2). The
    # paper reports a single shared sigma_add and sigma_prop for both
    # analytes; in nlmixr2 each endpoint needs its own residual parameter,
    # so the parent (Cc) uses the bare-named propSd / addSd and the
    # metabolite (Cc_9oh) uses propSd_9oh / addSd_9oh -- both endpoints
    # are initialised with the SAME Table 2 / A2 estimates to faithfully
    # encode the paper's shared-residual structure.
    propSd     <- 0.297; label("Proportional residual error on risperidone plasma concentration (fraction)") # Table A2 (Phase-3): sigma_prop = 29.7%
    addSd      <- 0.137; label("Additive residual error on risperidone plasma concentration (ng/mL)")        # Table A2 (Phase-3): sigma_add  = 0.137 ng/mL
    propSd_9oh <- 0.297; label("Proportional residual error on 9-OH-risperidone plasma concentration (fraction; shared with risperidone)") # Table A2: same sigma_prop = 29.7%
    addSd_9oh  <- 0.137; label("Additive residual error on 9-OH-risperidone plasma concentration (ng/mL; shared with risperidone)")        # Table A2: same sigma_add  = 0.137 ng/mL

    # ------------------- PANSS PD structural parameters ---------------------
    # Linear-scale parameters with additive normally-distributed IIV (Table
    # 3 footnote a: 'All are additive variability parameters, hence
    # reporting standard deviation in the same unit as parameter'). TPROG
    # is in weeks; the model() block converts model time (in hours from
    # first SC dose) to weeks via T_weeks = t / 168.
    bsl   <- 94.9;  label("Baseline PANSS total score (PANSS units)")                                          # Table 3 exposure-response: BSL = 94.9 (95% CI 93.9-95.8)
    pmax  <- 0.06;  label("Maximum relative placebo decrease in PANSS from baseline (fraction)")                # Table 3 exposure-response: PMAX = 0.06 (95% CI 0.04-0.09)
    tprog <- 1.7;   label("Time to 63.2% of the maximum placebo effect (weeks)")                                # Table 3 exposure-response: TPROG = 1.7 weeks (95% CI 1.06-2.0)
    pow   <- 2.1;   label("Weibull shape parameter on the placebo time course (unitless)")                      # Table 3 exposure-response: POW = 2.1 (95% CI 1.1-10.1)
    drift <- -1.2;  label("Linear drift in PANSS over time (PANSS units per week)")                             # Table 3 exposure-response: DRIFT = -1.2 (95% CI -1.5, -1.0)
    emax  <- 0.054; label("Maximum drug-effect relative decrease in PANSS from baseline (fraction)")            # Table 3 exposure-response: Emax = 0.054 (95% CI 0.02-0.09)
    ec50  <- 4.6;   label("Total active moiety concentration eliciting half-maximum drug effect (ng/mL)")       # Table 3 exposure-response: EC50 = 4.6 ng/mL (95% CI 0.18-33.4)

    # ------------------- IIV on PANSS PD ------------------------------------
    # Additive normal etas; omega = SD reported directly in Table 3 with the
    # same units as the corresponding parameter. No IIV reported for TPROG,
    # POW, Emax, or EC50 ('the IIV in EC50 could not be estimated, probably
    # due to the sparseness in this design, and was set to zero').
    etabsl   ~ 49.00     # Table 3: BSV(BSL)   = 7.0 (SD, PANSS units); variance = 7.0^2  = 49.00
    etapmax  ~ 0.000196  # Table 3: BSV(PMAX)  = 0.014 (SD, fraction);   variance = 0.014^2 = 0.000196
    etadrift ~ 1.69      # Table 3: BSV(DRIFT) = 1.3 (SD, PANSS units / week); variance = 1.3^2 = 1.69

    # ------------------- Residual error on PANSS ----------------------------
    addSd_PANSS <- 5.5; label("Additive residual error on PANSS total score (PANSS units)") # Table 3 exposure-response: RUV = 5.5 PANSS units
  })

  model({
    # Time-conversion. The PK rate constants are in 1/h and the dosing
    # interval is 28 days = 672 h; the PANSS PD equation uses time in
    # weeks (TPROG and DRIFT are reported per week). t is the rxode2
    # model time relative to the first dose, in hours; t_weeks is the
    # corresponding time in weeks.
    t_weeks <- t / 168

    # ------------------- Individual PK parameters ---------------------------
    ka1     <- exp(lka1     + etalka1)
    ka2     <- exp(lka2     + etalka2)
    ktr     <- exp(lktr     + etalktr)
    kel     <- exp(lkel)                       # no IIV reported on krel
    kel_9oh <- exp(lkel_9oh + etalkel_9oh)
    k12     <- exp(lk12     + etalk12)
    k21     <- exp(lk21     + etalk21)
    vc      <- exp(lvc      + etalvc)

    # CYP2D6 phenotype shift on the metabolite formation rate kr9.
    # Pooled EM + Inconclusive reference: both CYP2D6_PM = 0 and
    # CYP2D6_IM = 0 -> phenotype_factor = 1; PM -> 1 + (-0.94) = 0.06;
    # IM -> 1 + (-0.76) = 0.24.
    cyp2d6_factor <- 1 + e_cyp2d6_im_kmet * CYP2D6_IM + e_cyp2d6_pm_kmet * CYP2D6_PM
    kmet          <- exp(lkmet + etalkmet) * cyp2d6_factor

    # ------------------- Individual PANSS PD parameters ---------------------
    bsl_i   <- bsl   + etabsl
    pmax_i  <- pmax  + etapmax
    drift_i <- drift + etadrift

    # ------------------- ODE system -----------------------------------------
    # Compartments (in declaration order so the rxode2 cmt slots are stable
    # against any explicit cmt() in the event table):
    #   1. depot       - SC injection site (receives full dose); single
    #                    depot encoding with two parallel exit rates:
    #                    ka1 directly to central (fast peak path) and
    #                    ktr into the transit chain (slow ATRIGEL release).
    #                    This is the simplest mass-conserving interpretation
    #                    of Figure 1 of the source paper that uses only the
    #                    rate constants reported in Table 2; the implicit
    #                    fast-fraction is ka1 / (ka1 + ktr) ~ 18% and the
    #                    transit-fraction is ktr / (ka1 + ktr) ~ 82%. The
    #                    paper inherits the structural model from Gomeni 2013
    #                    and Laffont 2014/2015 without reporting an explicit
    #                    bioavailability split; see vignette Assumptions and
    #                    deviations.
    #   2-6. transit1..transit5 - 5-compartment transit chain; the slow
    #                    ATRIGEL release path enters transit1 at rate ktr
    #                    from depot, transits through the chain at rate ktr
    #                    per stage, and is absorbed from transit5 to central
    #                    at rate ka2.
    #   7. central     - risperidone central compartment.
    #   8. peripheral1 - risperidone peripheral compartment (k12 / k21
    #                    exchange).
    #   9. central_9oh - 9-OH-risperidone central compartment (apparent
    #                    volume of distribution constrained equal to V per
    #                    the paper's identifiability note).
    d/dt(depot)       <- -(ka1 + ktr) * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ka2 * transit5
    d/dt(central)     <-  ka1 * depot + ka2 * transit5 - (kel + kmet + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_9oh) <-  kmet * central - kel_9oh * central_9oh

    # ------------------- Observation variables ------------------------------
    # Plasma concentrations in ng/mL. Dose in mg, volumes in L: amount / V
    # gives mg/L = ug/mL; multiply by 1000 to obtain ng/mL, matching the
    # Table 1 / Figure 2 observed concentration units.
    Cc     <- 1000 * central     / vc
    Cc_9oh <- 1000 * central_9oh / vc

    # Total active moiety in risperidone-equivalent ng/mL. The molecular-
    # weight correction 410/426 converts 9-OH-risperidone mass concentration
    # to risperidone-equivalent concentration per the paper's PANSS PK/PD
    # methodology: '[AM] = [risperidone] + [9-OH-risperidone] * (410/426)'.
    # The rounded ratios 410 (risperidone MW = 410.49 g/mol) and 426 (9-OH-
    # risperidone MW = 426.49 g/mol) are used verbatim as in the source.
    AM <- Cc + Cc_9oh * (410 / 426)

    # PANSS PK/PD equation (Figure 1 box and Results 'PK/PD analysis for
    # PANSS score'): Predicted PANSS = BSL * (1 - PMAX * (1 - exp(-(T /
    # TPROG)^POW)) - DrugEff) + DRIFT * T, where DrugEff = Emax * AM /
    # (EC50 + AM). T is in weeks.
    DrugEff <- emax * AM / (ec50 + AM)
    PANSS   <- bsl_i * (1 - pmax_i * (1 - exp(-(t_weeks / tprog)^pow)) - DrugEff) + drift_i * t_weeks

    # ------------------- Residual error -------------------------------------
    # Plasma concentrations: combined additive + proportional error common
    # to risperidone and 9-OH-risperidone (Methods 'Pharmacokinetic model'
    # and Table 2 sigma_add / sigma_prop). PANSS: additive error.
    Cc     ~ add(addSd)     + prop(propSd)
    Cc_9oh ~ add(addSd_9oh) + prop(propSd_9oh)
    PANSS  ~ add(addSd_PANSS)
  })
}
