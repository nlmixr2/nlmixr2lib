Zhou_2018_remimazolam_moaas <- function() {
  description <- paste(
    "Sequential population PK/PD ordered-categorical (proportional-odds)",
    "model for the Modified Observer's Assessment of Alertness/Sedation",
    "(MOAA/S) score produced by remimazolam tosilate (development code",
    "HR7056) in 63 Chinese healthy adult volunteers after a single 1-minute",
    "intravenous injection of 0.01 to 0.45 mg/kg (Zhou 2018, NCT01970072).",
    "The three-compartment arterial-plasma pharmacokinetic model of",
    "modellib('Zhou_2018_remimazolam') is carried forward unchanged and",
    "drives a massless effect compartment with its own equilibration rate",
    "constant ke0, slower than the one fit to the Bispectral Index in the",
    "companion modellib('Zhou_2018_remimazolam_bis'). The drug raises each",
    "of five baseline cumulative logits by a common sigmoid Emax term in",
    "effect-site concentration; inverse-logit transformation gives the",
    "cumulative probabilities P(MOAA/S <= x) for x = 0 to 4, with",
    "P(MOAA/S <= 5) fixed at 1, and the six category probabilities follow by",
    "subtraction. The primary observed output is P(MOAA/S <= 1), the paper's",
    "own target effect for loss of consciousness. Age and sex were screened",
    "on the pharmacodynamic parameters and neither was retained.",
    sep = " "
  )
  reference <- paste(
    "Zhou Y, Hu P, Huang Y, Sang N, Song K, Wang H, Wen J, Jiang J, Chen X.",
    "Population Pharmacokinetic/Pharmacodynamic Model-Guided Dosing",
    "Optimization of a Novel Sedative HR7056 in Chinese Healthy Subjects.",
    "Front Pharmacol. 2018;9:1316. doi:10.3389/fphar.2018.01316. PMC6252322.",
    "Pharmacokinetic layer: see modellib('Zhou_2018_remimazolam').",
    sep = " "
  )
  vignette <- "Zhou_2018_remimazolam"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "ng/mL (remimazolam; the model outputs are MOAA/S category probabilities in 0-1)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhou 2018 Methods: arterial plasma
  # sampling for the PK states, and the 'link' model equation
  # dCe/dt = keo * (Cp - Ce), which makes the effect state a concentration
  # rather than an amount. The category probabilities are algebraic in the
  # effect-site concentration and carry no ODE state of their own.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "remimazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  # Screened on the MOAA/S pharmacodynamic parameters but NOT retained:
  # 'Covariate analysis of the effects of gender and age on parameters had
  # been undertaken for HR7056. However, no covariate effects were
  # statistically significant at the 0.01 level' (Zhou 2018 Results,
  # Population Pharmacodynamic Modeling For MOAA/S). Weight and height were
  # screened on the PK and BIS models but the MOAA/S section names only age
  # and gender, so only those two are recorded here. No point estimates are
  # reported for either.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      notes = "Screened on the MOAA/S pharmacodynamic parameters; not retained. Cohort median 27 years (range 18-44)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the MOAA/S pharmacodynamic parameters; not retained. The paper reports the covariate as 'gender' without stating the coding; SEXF is the canonical polarity (1 = female). 12 of 63 subjects were female."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    age_range = "18-44 years",
    age_median = "27 years",
    weight_range = "52.8-83.8 kg",
    weight_median = "63.8 kg",
    height_range = "156-184 cm",
    height_median = "169 cm",
    sex_female_pct = 100 * 12 / 63,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers. Eligibility: ages 18-55 years inclusive, weight 50-100 kg, BMI 18-26 kg/m^2.",
    dose_range = "Single 1-minute intravenous injection of 0.01, 0.02, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35 or 0.45 mg/kg (11 ascending-dose cohorts).",
    regions = "China (single centre; Clinical Pharmacology Research Center, Peking Union Medical College Hospital, Beijing)",
    n_observations = "1197 MOAA/S scores, recorded pre-dose and at 1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60 and 120 min after the start of the injection, alongside 1197 arterial plasma concentrations.",
    biomarkers = "Modified Observer's Assessment of Alertness/Sedation (MOAA/S), a six-level ordinal sedation scale scored 5 (fully alert, responds readily to name spoken in a normal tone) down to 0 (does not respond to a painful trapezius squeeze); Zhou 2018 Table S1.",
    notes = "Single-centre, double-blinded, randomised, single-ascending-dose study registered as NCT01970072. PK and PD were fit in two stages: individual PK parameters were estimated first with the three-compartment model of Table 1 and used to simulate the arterial concentration profiles that drive this PD model (Discussion). The categorical MOAA/S analysis was run in NONMEM 7.2 with the conditional Laplacian method, unlike the BIS analysis which used Phoenix NLME 1.2. The study's stopping rule was loss of consciousness (MOAA/S < 2) for at least 5 min in more than 50% of a cohort. Estimates from Zhou 2018 Table 2, 'FINAL MOAA/S MODEL' block."
  )

  ini({
    # ---------------------------------------------------------------------
    # Pharmacokinetic layer, Zhou 2018 Table 1. Identical to
    # modellib('Zhou_2018_remimazolam'); reproduced here so this PD model is
    # self-contained. The two-stage fit (Discussion) means these values were
    # not re-estimated during the PD step, but the paper does not mark them
    # FIX either -- they are the converged Table 1 estimates carried forward,
    # so they are encoded as ordinary estimated parameters exactly as in the
    # PK file.
    # ---------------------------------------------------------------------
    lvc <- log(2.11)
    label("Central (arterial) volume of distribution, V1 (L)") # Zhou 2018 Table 1: Central (arterial) volume = 2.11 L (RSE 4.0%)
    lcl <- log(1.49)
    label("Elimination clearance, CL (L/min)") # Zhou 2018 Table 1: Elimination clearance = 1.49 L/min (RSE 1.9%)
    lvp <- log(10.5)
    label("First peripheral volume of distribution, V2 (L)") # Zhou 2018 Table 1: Peripheral volume V2 = 10.5 L (RSE 3.9%)
    lq <- log(0.96)
    label("Inter-compartmental clearance to the first peripheral compartment, Cl2 (L/min)") # Zhou 2018 Table 1: Inter-compartmental clearance Cl2 = 0.96 L/min (RSE 3.7%)
    lvp2 <- log(22.7)
    label("Second peripheral volume of distribution, V3 (L)") # Zhou 2018 Table 1: Peripheral volume V3 = 22.7 L (RSE 4.1%)
    lq2 <- log(0.27)
    label("Inter-compartmental clearance to the second peripheral compartment, Cl3 (L/min)") # Zhou 2018 Table 1: Inter-compartmental clearance Cl3 = 0.27 L/min (RSE 4.0%)

    # ---------------------------------------------------------------------
    # Baseline cumulative logits, Zhou 2018 Table 2 rows B1 to B5
    # ('B1 ~ B5, Baseline value'). The Methods define
    #   Logit(x) = Baseline(x) + IMax * Ce^gamma / (IC50^gamma + Ce^gamma) + eta
    #   P(x)     = 1 / (1 + exp(-Logit(x)))
    # as the cumulative probability 'that the MOAA/S score was less than or
    # equal to that category', with 'the probability to observe a MOAA/S
    # score <= 5 being 1'. Six categories (0-5) with the top cumulative
    # probability fixed at 1 leave exactly five free cutpoints, x = 0 to 4,
    # which is the number of baselines reported. The tabulated values
    # increase monotonically in table order (-8.52 < -2.44 < -1.30 < -1.15 <
    # -0.895), as cumulative logits must, which is what fixes B1 to the
    # lowest category and B5 to the highest.
    # ---------------------------------------------------------------------
    b1 <- -8.52
    label("Baseline cumulative logit for MOAA/S <= 0 (logit scale)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: B1 = -8.52 (RSE 13.8%)
    b2 <- -2.44
    label("Baseline cumulative logit for MOAA/S <= 1 (logit scale)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: B2 = -2.44 (RSE 10.5%)
    b3 <- -1.30
    label("Baseline cumulative logit for MOAA/S <= 2 (logit scale)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: B3 = -1.30 (RSE 13.0%)
    b4 <- -1.15
    label("Baseline cumulative logit for MOAA/S <= 3 (logit scale)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: B4 = -1.15 (RSE 16.8%)
    b5 <- -0.895
    label("Baseline cumulative logit for MOAA/S <= 4 (logit scale)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: B5 = -0.895 (RSE 21.6%)

    # ---------------------------------------------------------------------
    # Drug effect on the cumulative logits, Zhou 2018 Table 2 ('FINAL MOAA/S
    # MODEL' block). The same sigmoid Emax term is added to every cumulative
    # logit, i.e. a proportional-odds structure: 'The drug acts to increase
    # the baseline logit according to a sigmoid Emax model based on the
    # predicted effect-site concentration (Ce)' (Results). Imax is therefore
    # a logit increment, not a probability, despite the Methods glossing it
    # as 'the maximal achievable probability'.
    #
    # The effect compartment has its OWN equilibration rate constant here
    # (0.05 1/min) even though the Methods state that 'MOAA/S score related
    # effect-site concentration was the same as that of BIS' -- Table 2
    # reports 0.0855 1/min for BIS and 0.05 1/min for MOAA/S, and the
    # Discussion quotes the two corresponding half-lives separately (8.1 min
    # and 13.9 min). The tabulated values are used.
    # ---------------------------------------------------------------------
    lke0 <- log(0.05)
    label("Effect-compartment equilibration rate constant, Ke0 (1/min)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: Ke0 = 0.05 1/min (RSE 7.3%); Discussion reports the corresponding half-life as 13.9 min (log(2)/0.05 = 13.86)
    lec50 <- log(436)
    label("Effect-site concentration producing half of Imax, IC50 (ng/mL)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: IC50 = 436 ng/mL (RSE 15.6%)
    lhill <- log(1.50)
    label("Hill coefficient of the sigmoid Emax term on the cumulative logits (unitless)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: Hill coefficient gamma = 1.50 (RSE 14.6%)
    limax <- log(27.9)
    label("Maximum drug-induced increase in each cumulative logit, Imax (logit units)") # Zhou 2018 Table 2 FINAL MOAA/S MODEL: Imax = 27.9 (RSE 9.4%)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 2's 'IIV (RSE%)' column follows the
    # same convention the paper states for the PK model: theta_i = theta_TV *
    # exp(eta_i) with the tabulated figure being omega * 100, the SD of eta
    # in the log domain. Internal variance is (IIV% / 100)^2. Table 2 shows
    # '-' (not applicable) for B1 to B5, so the baselines carry no IIV; the
    # additive 'eta' written into the Methods logit equation is taken to be
    # realised through the log-normal IIV on the four drug-effect parameters
    # that Table 2 does quantify (see the vignette 'Assumptions and
    # deviations').
    #
    # PK-layer IIV is reproduced from Table 1 so the model can be simulated
    # end-to-end.
    # ---------------------------------------------------------------------
    etalvc ~ 0.0196 # Zhou 2018 Table 1: IIV on central volume = 14.0% (RSE 5.6%); omega^2 = 0.140^2
    etalcl ~ 0.013225 # Zhou 2018 Table 1: IIV on elimination clearance = 11.5% (RSE 2.3%); omega^2 = 0.115^2
    etalvp ~ 0.014884 # Zhou 2018 Table 1: IIV on V2 = 12.2% (RSE 5.2%); omega^2 = 0.122^2
    etalq ~ 0.017689 # Zhou 2018 Table 1: IIV on Cl2 = 13.3% (RSE 4.5%); omega^2 = 0.133^2
    etalvp2 ~ 0.064009 # Zhou 2018 Table 1: IIV on V3 = 25.3% (RSE 5.8%); omega^2 = 0.253^2
    etalq2 ~ 0.034969 # Zhou 2018 Table 1: IIV on Cl3 = 18.7% (RSE 7.8%); omega^2 = 0.187^2

    etalke0 ~ 0.146689 # Zhou 2018 Table 2 FINAL MOAA/S MODEL: IIV on Ke0 = 38.3% (RSE 11.0%); omega^2 = 0.383^2
    etalec50 ~ 0.207025 # Zhou 2018 Table 2 FINAL MOAA/S MODEL: IIV on IC50 = 45.5% (RSE 29.8%); omega^2 = 0.455^2
    etalhill ~ 0.309136 # Zhou 2018 Table 2 FINAL MOAA/S MODEL: IIV on the Hill coefficient = 55.6% (RSE 40.1%); omega^2 = 0.556^2
    etalimax ~ 0.030976 # Zhou 2018 Table 2 FINAL MOAA/S MODEL: IIV on Imax = 17.6% (RSE 32.6%); omega^2 = 0.176^2

    # ---------------------------------------------------------------------
    # Placeholder residual error. Zhou 2018 fit the MOAA/S data by an exact
    # ordered-categorical likelihood (NONMEM 7.2, conditional Laplacian) and
    # therefore estimates no residual error; Table 2 prints a sigma for the
    # BIS model only. rxode2 / nlmixr2 require an observation model, so the
    # typical-value probability is exposed with a small placeholder additive
    # residual, following the Koo_2012_propofol.R and Shin_2014_sevoflurane.R
    # precedent for Bernoulli-likelihood probability outputs. This value is
    # NOT from the source paper.
    # ---------------------------------------------------------------------
    addSd_prob_moaas_le1 <- 0.05
    label("Additive residual error on the probability of MOAA/S <= 1 (unitless; placeholder)")
  })

  model({
    # Amounts are carried in mg and volumes in L, so central/vc is mg/L;
    # Zhou 2018 reports IC50 and the assay range in ng/mL, so concentrations
    # are converted with 1 mg/L = 1000 ng/mL.
    mgL_to_ngmL <- 1000

    # Individual disposition parameters.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Individual pharmacodynamic parameters.
    ke0 <- exp(lke0 + etalke0)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)
    imax <- exp(limax + etalimax)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition with intravenous input into the central
    # compartment (Zhou 2018 Table 1).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Arterial plasma concentration (ng/mL).
    Cc <- central / vc * mgL_to_ngmL

    # Massless effect compartment, Zhou 2018 Methods: dCe/dt = keo * (Cp -
    # Ce). The state is a concentration in ng/mL, not an amount.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Common sigmoid Emax increment added to every cumulative logit
    # (proportional-odds structure).
    drug <- imax * effect^hill / (ec50^hill + effect^hill)

    # Cumulative logits. The baselines b2 to b5 are wrapped in parentheses so
    # rxode2's mu-reference parser does not read 'b<n> + drug' as a
    # population-parameter-plus-eta pair.
    logit_moaas_le0 <- b1 + drug
    logit_moaas_le1 <- (b2) + drug
    logit_moaas_le2 <- (b3) + drug
    logit_moaas_le3 <- (b4) + drug
    logit_moaas_le4 <- (b5) + drug

    # Cumulative probabilities P(MOAA/S <= x), Zhou 2018 Methods:
    # P(x) = 1 / (1 + exp(-Logit(x))). P(MOAA/S <= 5) is fixed at 1 and so
    # has no cutpoint of its own.
    prob_moaas_le0 <- expit(logit_moaas_le0)
    prob_moaas_le1 <- expit(logit_moaas_le1)
    prob_moaas_le2 <- expit(logit_moaas_le2)
    prob_moaas_le3 <- expit(logit_moaas_le3)
    prob_moaas_le4 <- expit(logit_moaas_le4)

    # Category probabilities by subtraction from the cumulative
    # probabilities, Zhou 2018 Methods.
    prob_moaas0 <- prob_moaas_le0
    prob_moaas1 <- prob_moaas_le1 - prob_moaas_le0
    prob_moaas2 <- prob_moaas_le2 - prob_moaas_le1
    prob_moaas3 <- prob_moaas_le3 - prob_moaas_le2
    prob_moaas4 <- prob_moaas_le4 - prob_moaas_le3
    prob_moaas5 <- 1 - prob_moaas_le4

    # Expected MOAA/S score, a convenience summary of the six category
    # probabilities. Not a quantity the source paper reports.
    moaas_expected <- 0 * prob_moaas0 + 1 * prob_moaas1 + 2 * prob_moaas2 +
      3 * prob_moaas3 + 4 * prob_moaas4 + 5 * prob_moaas5

    # Observed output: the probability of MOAA/S < 2, which is
    # P(MOAA/S <= 1) exactly. This is the paper's own target effect for loss
    # of consciousness -- 'MOAA/S score < 2 for at least 1.5 h' -- and the
    # P0_1 trace of Figure 5B. Placeholder residual; see the ini() comment.
    prob_moaas_le1 ~ add(addSd_prob_moaas_le1)
  })
}
