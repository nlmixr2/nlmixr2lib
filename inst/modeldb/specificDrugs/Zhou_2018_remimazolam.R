Zhou_2018_remimazolam <- function() {
  description <- paste(
    "Three-compartment population pharmacokinetic model for remimazolam",
    "tosilate (development code HR7056) in 63 Chinese healthy adult",
    "volunteers given a single 1-minute intravenous injection of 0.01 to",
    "0.45 mg/kg across 11 ascending-dose cohorts (Zhou 2018, NCT01970072).",
    "Disposition is parameterised as clearances and volumes, with linear",
    "elimination from the central compartment and two peripheral",
    "compartments; the central volume is an arterial central volume because",
    "the study sampled arterial rather than venous plasma. Body weight,",
    "height, age and sex were each screened as covariates on every",
    "pharmacokinetic parameter and none reached the 0.01 significance level,",
    "so the final model carries no covariates. Inter-individual variability",
    "is log-normal on all six disposition parameters and the residual error",
    "is proportional. This is the pharmacokinetic layer of a two-stage",
    "PK/PD analysis; the sedation models fit against it are",
    "modellib('Zhou_2018_remimazolam_bis') and",
    "modellib('Zhou_2018_remimazolam_moaas').",
    sep = " "
  )
  reference <- paste(
    "Zhou Y, Hu P, Huang Y, Sang N, Song K, Wang H, Wen J, Jiang J, Chen X.",
    "Population Pharmacokinetic/Pharmacodynamic Model-Guided Dosing",
    "Optimization of a Novel Sedative HR7056 in Chinese Healthy Subjects.",
    "Front Pharmacol. 2018;9:1316. doi:10.3389/fphar.2018.01316. PMC6252322.",
    sep = " "
  )
  vignette <- "Zhou_2018_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhou 2018 Methods ('Arterial plasma
  # samples were collected...') and the Table 1 parameterisation, which names
  # a central (arterial) volume plus peripheral volumes V2 and V3.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Screened in turn on every pharmacokinetic parameter but NOT retained:
  # 'Covariate analysis of the possible effects of body weight, height, age,
  # and gender was examined for each pharmacokinetic parameter. No covariate
  # effects were statistically significant at the 0.01 level (dOFV = -6.63)'
  # (Zhou 2018 Results, Population Pharmacokinetic Analysis of HR7056). No
  # point estimates are reported for any of them, so they are documented here
  # rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Screened on every PK parameter; not retained. Cohort median 63.8 kg (range 52.8-83.8); eligibility required 50-100 kg with BMI 18-26 kg/m^2."
    ),
    HT = list(
      description = "Body height at baseline.",
      units = "cm",
      type = "continuous",
      notes = "Screened on every PK parameter; not retained. Cohort median 169 cm (range 156-184)."
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      notes = "Screened on every PK parameter; not retained. Cohort median 27 years (range 18-44); eligibility 18-55 years."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on every PK parameter; not retained. The paper reports the covariate as 'gender' without stating the coding; SEXF is the canonical polarity (1 = female). 12 of 63 subjects were female. The authors attribute the absence of covariate effects to the strict entry criteria and the resulting narrow demographic spread (Discussion)."
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
    disease_state = "Healthy volunteers. Eligibility: ages 18-55 years inclusive, weight 50-100 kg, BMI 18-26 kg/m^2. 194 screened subjects failed the inclusion criteria, most commonly for abnormal laboratory or physical findings, a positive smoking test, or an ineligible age or BMI.",
    dose_range = "Single 1-minute intravenous injection of 0.01, 0.02, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35 or 0.45 mg/kg (11 ascending-dose cohorts).",
    regions = "China (single centre; Clinical Pharmacology Research Center, Peking Union Medical College Hospital, Beijing)",
    n_observations = "1197 arterial plasma concentrations (and, in the companion PD models, 1197 BIS values and 1197 MOAA/S scores).",
    notes = "Single-centre, double-blinded, randomised, single-ascending-dose study registered as NCT01970072. A further 16 subjects received midazolam (6 at 0.075 mg/kg, 10 at 0.12 mg/kg) as an active comparator and were not part of this model. Arterial plasma was sampled pre-dose and at 1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 45, 60, 120, 180, 240, 360 and 480 min post-dose and assayed by UPLC-MS/MS validated over 0.5-1000 ng/mL (LLOQ 0.5 ng/mL). Fitted in Phoenix NLME 1.2 with FOCE-I. Demographics from Zhou 2018 Results; parameter estimates from Table 1. Total volume of distribution 2.11 + 10.5 + 22.7 = 35.3 L, which the Discussion compares with 34.9 L for remimazolam in Wiltshire 2012."
  )

  ini({
    # ---------------------------------------------------------------------
    # Disposition typical values, Zhou 2018 Table 1 ('Typical value (RSE%)'
    # column). The paper parameterises the three-compartment model as
    # clearances and volumes; V2/Cl2 map to the canonical peripheral1 pair
    # and V3/Cl3 to the peripheral2 pair. 95% CIs in Table 1 are from
    # log-likelihood profiling (Statistical Analysis section).
    # ---------------------------------------------------------------------
    lvc <- log(2.11)
    label("Central (arterial) volume of distribution, V1 (L)") # Zhou 2018 Table 1: Central (arterial) volume = 2.11 L (RSE 4.0%), 95% CI 1.94-2.27
    lcl <- log(1.49)
    label("Elimination clearance, CL (L/min)") # Zhou 2018 Table 1: Elimination clearance = 1.49 L/min (RSE 1.9%), 95% CI 1.44-1.55
    lvp <- log(10.5)
    label("First peripheral volume of distribution, V2 (L)") # Zhou 2018 Table 1: Peripheral volume V2 = 10.5 L (RSE 3.9%), 95% CI 9.73-11.3
    lq <- log(0.96)
    label("Inter-compartmental clearance to the first peripheral compartment, Cl2 (L/min)") # Zhou 2018 Table 1: Inter-compartmental clearance Cl2 = 0.96 L/min (RSE 3.7%), 95% CI 0.89-1.03
    lvp2 <- log(22.7)
    label("Second peripheral volume of distribution, V3 (L)") # Zhou 2018 Table 1: Peripheral volume V3 = 22.7 L (RSE 4.1%), 95% CI 20.9-24.6
    lq2 <- log(0.27)
    label("Inter-compartmental clearance to the second peripheral compartment, Cl3 (L/min)") # Zhou 2018 Table 1: Inter-compartmental clearance Cl3 = 0.27 L/min (RSE 4.0%), 95% CI 0.250-0.295

    # ---------------------------------------------------------------------
    # Inter-individual variability, Zhou 2018 Table 1 ('IIV% (RSE%)' column).
    # The paper states the convention explicitly in the Population
    # Pharmacokinetic Modeling section: theta_i = theta_TV * exp(eta_i), and
    # 'Inter-individual variability is reported as omega, the SD of eta in
    # the log domain, which is approximately the coefficient of variation in
    # the standard domain'. The tabulated IIV% is therefore omega * 100 and
    # the internal variance is (IIV% / 100)^2 -- NOT log(CV^2 + 1), which
    # would double-count the log-normal transformation the paper has already
    # applied.
    # ---------------------------------------------------------------------
    etalvc ~ 0.0196 # Zhou 2018 Table 1: IIV on central volume = 14.0% (RSE 5.6%); omega^2 = 0.140^2
    etalcl ~ 0.013225 # Zhou 2018 Table 1: IIV on elimination clearance = 11.5% (RSE 2.3%); omega^2 = 0.115^2
    etalvp ~ 0.014884 # Zhou 2018 Table 1: IIV on V2 = 12.2% (RSE 5.2%); omega^2 = 0.122^2
    etalq ~ 0.017689 # Zhou 2018 Table 1: IIV on Cl2 = 13.3% (RSE 4.5%); omega^2 = 0.133^2
    etalvp2 ~ 0.064009 # Zhou 2018 Table 1: IIV on V3 = 25.3% (RSE 5.8%); omega^2 = 0.253^2
    etalq2 ~ 0.034969 # Zhou 2018 Table 1: IIV on Cl3 = 18.7% (RSE 7.8%); omega^2 = 0.187^2

    # ---------------------------------------------------------------------
    # Residual error. 'A proportional error model was used for the residual
    # random effects' (Methods), and the Results restate Table 1's sigma as a
    # percentage: 'residual variability was 13.8%'. That restatement is what
    # pins sigma to the proportional SD rather than to a variance.
    # ---------------------------------------------------------------------
    propSd <- 0.138
    label("Proportional residual error on remimazolam arterial plasma concentration (fraction)") # Zhou 2018 Table 1: sigma = 0.138 (RSE 2.9%), 95% CI 0.130-0.145; Results: 'residual variability was 13.8%'
  })

  model({
    # Amounts are carried in mg and volumes in L, so central/vc is mg/L;
    # Zhou 2018 reports every concentration (the 0.5-1000 ng/mL assay range,
    # and the IC50 values of the companion PD models) in ng/mL, so
    # concentrations are converted with 1 mg/L = 1000 ng/mL.
    mgL_to_ngmL <- 1000

    # Individual disposition parameters.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition with intravenous input into the central
    # compartment (Zhou 2018 Methods: HR7056 was 'administered as a 1-min IV
    # injection').
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Arterial plasma concentration (ng/mL).
    Cc <- central / vc * mgL_to_ngmL

    Cc ~ prop(propSd)
  })
}
