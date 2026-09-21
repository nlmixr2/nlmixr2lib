Zhou_2018_remimazolam_bis <- function() {
  description <- paste(
    "Sequential population PK/PD model for the Bispectral Index (BIS) depth",
    "of sedation produced by remimazolam tosilate (development code HR7056)",
    "in 63 Chinese healthy adult volunteers after a single 1-minute",
    "intravenous injection of 0.01 to 0.45 mg/kg (Zhou 2018, NCT01970072).",
    "The three-compartment arterial-plasma pharmacokinetic model of",
    "modellib('Zhou_2018_remimazolam') is carried forward unchanged (the",
    "authors fit PK and PD in two stages) and drives a massless effect",
    "compartment through the equilibration rate constant ke0. Effect-site",
    "concentration then enters an inhibitory sigmoid Imax model on BIS,",
    "E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma). Age, weight, height",
    "and sex were screened on the pharmacodynamic parameters and none was",
    "retained. Inter-individual variability is log-normal on ke0, IC50, the",
    "Hill coefficient, E0 and Imax. The companion ordered-categorical",
    "sedation-score model is modellib('Zhou_2018_remimazolam_moaas').",
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
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhou 2018 Methods: arterial plasma
  # sampling for the PK states, and the 'link' model equation
  # dCe/dt = keo * (Cp - Ce), which makes the effect state a concentration
  # rather than an amount.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "remimazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  # Screened on the pharmacodynamic parameters but NOT retained: 'inclusion of
  # age, weight, height, or gender did not result in a significant decrease in
  # the OFV. Therefore, no covariates were included in the final model' (Zhou
  # 2018 Results, Population Pharmacodynamic Modeling For BIS). The same four
  # were screened and rejected on the PK parameters; no point estimates are
  # reported for any of them.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Screened on the BIS pharmacodynamic parameters; not retained. Cohort median 63.8 kg (range 52.8-83.8)."
    ),
    HT = list(
      description = "Body height at baseline.",
      units = "cm",
      type = "continuous",
      notes = "Screened on the BIS pharmacodynamic parameters; not retained. Cohort median 169 cm (range 156-184)."
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      notes = "Screened on the BIS pharmacodynamic parameters; not retained. Cohort median 27 years (range 18-44)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the BIS pharmacodynamic parameters; not retained. The paper reports the covariate as 'gender' without stating the coding; SEXF is the canonical polarity (1 = female). 12 of 63 subjects were female."
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
    n_observations = "1197 BIS values, recorded pre-dose and at 1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60 and 120 min after the start of the injection, alongside 1197 arterial plasma concentrations.",
    biomarkers = "Bispectral Index (BIS), a processed-EEG index of sedation depth on a 0-100 scale where 100 indicates a fully awake subject.",
    notes = "Single-centre, double-blinded, randomised, single-ascending-dose study registered as NCT01970072. PK and PD were fit in two stages: individual PK parameters were estimated first with the three-compartment model of Table 1 and used to simulate the arterial concentration profiles that drive this PD model (Discussion). The BIS analysis was run in Phoenix NLME 1.2. Estimates from Zhou 2018 Table 2, 'FINAL BIS MODEL' block."
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
    # Pharmacodynamic parameters, Zhou 2018 Table 2 ('FINAL BIS MODEL'
    # block). The Methods give the structure as
    #   dCe/dt = keo * (Cp - Ce)   and   E = E0 - IMax * Ce^gamma /
    #                                        (IC50^gamma + Ce^gamma)
    # with E0 the drug-free BIS, IMax the maximum possible reduction in BIS,
    # IC50 the effect-site concentration halfway between E0 and E0 - IMax,
    # and gamma the Hill coefficient.
    # ---------------------------------------------------------------------
    lke0 <- log(0.0855)
    label("Effect-compartment equilibration rate constant, Ke0 (1/min)") # Zhou 2018 Table 2 FINAL BIS MODEL: Ke0 = 0.0855 1/min (RSE 7.4%); Discussion reports the corresponding half-life as 8.1 min (log(2)/0.0855 = 8.11)
    lec50 <- log(503)
    label("Effect-site concentration producing half of Imax, IC50 (ng/mL)") # Zhou 2018 Table 2 FINAL BIS MODEL: IC50 = 503 ng/mL (RSE 10.3%)
    lhill <- log(1.50)
    label("Hill coefficient of the inhibitory sigmoid Imax model (unitless)") # Zhou 2018 Table 2 FINAL BIS MODEL: Hill coefficient gamma = 1.50 (RSE 7.8%)
    lrbase <- log(95.3)
    label("Drug-free baseline Bispectral Index, E0 (BIS units)") # Zhou 2018 Table 2 FINAL BIS MODEL: E0 = 95.3 (RSE 0.4%)
    limax <- log(47.9)
    label("Maximum possible reduction in Bispectral Index, Imax (BIS units)") # Zhou 2018 Table 2 FINAL BIS MODEL: Imax = 47.9 (RSE 7.4%)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 2's 'IIV (RSE%)' column follows the
    # same convention the paper states for the PK model: theta_i = theta_TV *
    # exp(eta_i) with the tabulated figure being omega * 100, the SD of eta
    # in the log domain. The Simulations section confirms the distributional
    # family for the PD step -- 'the variability were randomly sampled from
    # the log-normal distributions obtained from the modeling'. Internal
    # variance is therefore (IIV% / 100)^2.
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

    etalke0 ~ 0.262144 # Zhou 2018 Table 2 FINAL BIS MODEL: IIV on Ke0 = 51.2% (RSE 25.1%); omega^2 = 0.512^2
    etalec50 ~ 0.168921 # Zhou 2018 Table 2 FINAL BIS MODEL: IIV on IC50 = 41.1% (RSE 20.2%); omega^2 = 0.411^2
    etalhill ~ 0.746496 # Zhou 2018 Table 2 FINAL BIS MODEL: IIV on the Hill coefficient = 86.4% (RSE 36.5%); omega^2 = 0.864^2
    etalrbase ~ 0.00024336 # Zhou 2018 Table 2 FINAL BIS MODEL: IIV on E0 = 1.56% (RSE 22.1%); omega^2 = 0.0156^2
    etalimax ~ 0.024336 # Zhou 2018 Table 2 FINAL BIS MODEL: IIV on Imax = 15.6% (RSE 30.1%); omega^2 = 0.156^2

    # ---------------------------------------------------------------------
    # Residual error on BIS. Zhou 2018 prints ONE residual parameter for the
    # BIS model, sigma = 0.0653, and the Methods describe the estimation as
    # using 'an additive residual error model'. Those two statements cannot
    # both be taken at face value: an additive SD of 0.0653 on a 0-100 BIS
    # scale would mean every BIS observation was predicted to within about a
    # tenth of a BIS unit, and it cannot produce the behaviour the paper's
    # own Discussion reports for its BIS visual predictive check -- 'the
    # simulated prediction presented a BIS higher than 100, which can be
    # ameliorated with a certain degree by changing the error model'. With
    # E0 = 95.3 and only 1.56% IIV on E0, reaching a simulated BIS above 100
    # requires several BIS units of residual spread at baseline, which is
    # what a 6.53% PROPORTIONAL residual delivers (0.0653 * 95.3 = 6.2 BIS
    # units) and what a 0.0653 additive residual cannot.
    #
    # Following the standing ambiguous-residual-error convention, both
    # readings are encoded, each carrying the single printed sigma: the
    # proportional term is the reading consistent with the paper's own
    # diagnostics and dominates everywhere on the BIS scale, and the additive
    # term is the paper's literal wording and acts only as a negligible floor
    # near BIS = 0. See the vignette 'Assumptions and deviations' section.
    # ---------------------------------------------------------------------
    propSd_BIS <- 0.0653
    label("Proportional residual error on BIS (fraction)") # Zhou 2018 Table 2 FINAL BIS MODEL: sigma = 0.0653 (RSE 2.2%), read as a proportional SD
    addSd_BIS <- 0.0653
    label("Additive residual error on BIS (BIS units)") # Zhou 2018 Table 2 FINAL BIS MODEL: sigma = 0.0653 (RSE 2.2%), read as an additive SD per the Methods wording
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
    rbase <- exp(lrbase + etalrbase)
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

    # Inhibitory sigmoid Imax model on BIS, Zhou 2018 Methods:
    # E = E0 - IMax * Ce^gamma / (IC50^gamma + Ce^gamma).
    BIS <- rbase - imax * effect^hill / (ec50^hill + effect^hill)

    # BIS is the only observed endpoint of this stage of the analysis: the
    # paper fit PK first and then fit the BIS data against the simulated
    # arterial concentration profiles (Discussion), so the plasma residual
    # error belongs to modellib('Zhou_2018_remimazolam') rather than here.
    # Cc remains available as a derived output column.
    BIS ~ prop(propSd_BIS) + add(addSd_BIS)
  })
}
