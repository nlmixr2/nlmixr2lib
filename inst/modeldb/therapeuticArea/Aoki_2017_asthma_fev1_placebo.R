Aoki_2017_asthma_fev1_placebo <- function() {
  description <- "Longitudinal placebo model of clinic pre-bronchodilator forced expiratory volume in 1 second (FEV1) in adults with inadequately controlled atopic asthma on budesonide/formoterol maintenance therapy, fit to the placebo and lowest-dose (AZD1981 10 mg twice daily) arms of the 12-week AZD1981 Phase IIb trial (NCT01197794; Aoki 2017). Baseline FEV1 is log-normal with linear-deviation effects of baseline FEV1 percent predicted and age and a multiplicative female-sex factor; the placebo effect is an additive step that is zero at the screening and baseline visits and constant at every post-randomization visit, with additive between-subject variability. Combined proportional and additive residual error. The model has no drug-effect term: the paper's candidate dose-effect functions are simulation devices with no AZD1981 estimates."
  reference <- paste(
    "Aoki Y, Roshammar D, Hamren B, Hooker AC. Model selection and averaging",
    "of nonlinear mixed-effect models for robust phase III dose selection.",
    "J Pharmacokinet Pharmacodyn. 2017;44(6):581-597.",
    "doi:10.1007/s10928-017-9550-0.",
    "Trial design and baseline demographics: Bateman ED, O'Brien C, Rugman P,",
    "Luke S, Ivanov S, Uddin M. Efficacy and safety of the CRTh2 antagonist",
    "AZD1981 as add-on therapy to inhaled corticosteroids and long-acting",
    "beta2-agonists in patients with atopic asthma. Drug Des Devel Ther.",
    "2018;12:1093-1106. doi:10.2147/DDDT.S147389.",
    sep = " "
  )
  vignette <- "Aoki_2017_asthma_fev1_placebo"
  units <- list(
    time = "week (time since the randomization / baseline visit, visit 2)",
    dosing = "n/a (placebo model; no dosing events)",
    concentration = "L (FEV1 absolute volume, observation FEV1)"
  )

  covariateData <- list(
    FEV1_PCTPRED = list(
      description = "Baseline (visit 2) FEV1 as percent of the predicted normal value",
      units = "% predicted",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear-deviation effect on baseline FEV1, (1 + e_fev1_pctpred_rbase * (FEV1_PCTPRED - 62.5)) (Aoki 2017 placebo-model equation; Table 1 theta3). Aoki 2017 defines the covariate as 'the percentage of FEV1 at visit 2 compared to the predicted normal' and centres it on the analysis-population mean, which is not printed in the paper or in the trial report (Bateman 2018). The reference of 62.5 % is the midpoint of the trial's 40-85 % predicted-normal inclusion range (Bateman 2018 Methods), a maintainer choice. The reference equation for predicted normal is not stated. Time-fixed.",
      source_name = "FEV1% of normal"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear-deviation effect on baseline FEV1, (1 + e_age_rbase * (AGE - 45)) (Aoki 2017 placebo-model equation; Table 1 theta4). Aoki 2017 centres age on the analysis-population mean, which is not printed; 45 years rounds the 45.4-year mean of the placebo and AZD1981 10 mg arms of the same trial (Bateman 2018 Table 1: 44.5 and 46.3 years), a maintainer choice. Time-fixed.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative factor on baseline FEV1: 1 for males and theta5 = 0.719 for females (Aoki 2017 placebo-model equation; Table 1 theta5), encoded as (1 - SEXF) + e_sexf_rbase * SEXF. The paper's indicator already has male as the reference, so no value inversion is needed.",
      source_name = "Sex (Male / Female)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 324L,
    n_studies = 1L,
    n_observations = 1803L,
    age_range = "adults >= 18 years (inclusion criterion); arm means 44.5 (placebo) and 46.3 (AZD1981 10 mg) years (Bateman 2018 Table 1)",
    sex_female_pct = 49.8,
    race_ethnicity = c(White = 72.0, Asian = 17.5, Black = 1.2, Other = 9.2),
    disease_state = "Inadequately controlled moderate-to-severe atopic asthma (prebronchodilator FEV1 40-85 % of predicted normal, >= 12 % reversibility, positive Phadiatop test, ACQ5 >= 1) on budesonide/formoterol 80/4.5 ug maintenance therapy",
    dose_range = "Placebo or AZD1981 10 mg twice daily orally for 12 weeks; the two arms were pooled as a 'placebo' dataset because they did not differ significantly",
    regions = "Europe, South and North America, Japan and South Africa (multicentre Phase IIb trial NCT01197794)",
    notes = "Aoki 2017 analysed the placebo and lowest-dose arms of the AZD1981 Phase IIb dose-finding trial (324 patients, 1803 clinic FEV1 measurements; seven visits, of which visit 1 is screening and visit 2 is baseline). Aoki 2017 does not tabulate demographics; age, sex and race are pooled from the placebo and 10 mg BID columns of Bateman 2018 Table 1, whose percentages imply n = 163 and n = 162 (325 in total, against the 324 Aoki analysed). Race: White 234, Asian 57, Black 4, other 30. Dosing records were not supplied to Aoki 2017."
  )

  ini({
    # Aoki 2017 Table 1 (NONMEM 7.3 estimates; placebo model fit to the
    # placebo and AZD1981 10 mg arms). The paper prints SDs for the random
    # effects; variances below are SD^2.
    lrbase <- log(2.51)
    label("Typical baseline FEV1 in a male patient at the reference age and FEV1 percent predicted (L)") # Table 1 theta2 'Baseline' 2.51 L (RSE 0.828%)
    plbmax <- 0.169
    label("Typical placebo effect on FEV1 at every post-baseline visit (L)") # Table 1 theta1 'Placebo effect' 0.169 L (RSE 11.9%)

    e_fev1_pctpred_rbase <- 0.0129
    label("Linear effect of baseline FEV1 percent predicted on baseline FEV1 (per percentage point)") # Table 1 theta3 'Covariate effect of FEV1% of normal' 0.0129 (RSE 3.87%); printed unit 'L-1'
    e_age_rbase <- -0.0105
    label("Linear effect of age on baseline FEV1 (1/year)") # Table 1 theta4 'Covariate effect of age' -0.0105 year-1 (RSE 5.01%)
    e_sexf_rbase <- 0.719
    label("Multiplicative factor on baseline FEV1 for female patients (unitless)") # Table 1 theta5 'Covariate effect of sex' 0.719 (RSE 1.43%)

    etaplbmax ~ 0.091809 # Table 1 omega1 'SD of IIV of the placebo effect' 0.303 L (RSE 8.88%); additive eta, variance 0.303^2
    etalrbase ~ 0.011025 # Table 1 omega2 'SD of IIV of the baseline' 0.105 (RSE 5.97%); exponential eta, variance 0.105^2

    propSd_FEV1 <- 0.0832
    label("Proportional residual error on FEV1 (fraction)") # Table 1 sigma1 'SD of proportional RUV' 0.0832 (RSE 7.60%)
    addSd_FEV1 <- 0.102
    label("Additive residual error on FEV1 (L)") # Table 1 sigma2 'SD of additive RUV' 0.102 L (RSE 18.8%)
  })

  model({
    # Baseline FEV1 (Aoki 2017 placebo-model equation):
    # theta2 * exp(eta2) * (1 + theta3 * (FEV1%normal - mean)) *
    # (1 + theta4 * (Age - mean)) * (1 if male, theta5 if female).
    # The population means are not printed; 62.5 % and 45 years are the
    # maintainers' reference values (see covariateData notes).
    rbase <- exp(lrbase + etalrbase) *
      (1 + e_fev1_pctpred_rbase * (FEV1_PCTPRED - 62.5)) *
      (1 + e_age_rbase * (AGE - 45)) *
      ((1 - SEXF) + e_sexf_rbase * SEXF)

    # Placebo effect: a step in visit, 0 at visits 1-2 (screening and
    # baseline) and theta1 + eta1 at visits 3-7. Time is weeks since the
    # baseline visit, so every post-baseline observation (t > 0) carries the
    # full placebo effect.
    plbi <- plbmax + etaplbmax
    plbon <- 0
    if (t > 0) {
      plbon <- 1
    }

    FEV1 <- rbase + plbi * plbon

    # NONMEM Y = IPRED * (1 + EPS1) + EPS2, i.e. combined proportional and
    # additive error with variance IPRED^2 * sigma1^2 + sigma2^2.
    FEV1 ~ add(addSd_FEV1) + prop(propSd_FEV1)
  })
}
