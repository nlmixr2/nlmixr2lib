You_2026_moxifloxacin_ttp14d <- function() {
  description <- "Biexponential time-to-positivity (TTP) early-bactericidal-activity model for the first 14 days of treatment in adults with newly diagnosed smear-positive pulmonary tuberculosis in the REMoxTB trial (control isoniazid + rifampicin + pyrazinamide + ethambutol, or moxifloxacin replacing ethambutol or isoniazid). MGIT liquid-culture TTP (h) rises with time on treatment as the sum of two exponentially growing terms; the second-term intercept and both slopes differ between the control and the two pooled moxifloxacin arms. Covariates: lung cavitation, sex and smoking on the first-term log amplitude; Asian or unreported race and sex on the second-term log amplitude; unreported race, age and smoking on the second slope. Exponential IIV on the first-term log amplitude and one shared IIV on both slopes; proportional residual error. A $PRED regression on time with no exposure term, no differential equations and no dosing events."
  reference <- paste(
    "You H, Simonsson USH. Applying a pharmacometrics-enabled machine",
    "learning analysis to predict 2-month culture conversion using phase 2a",
    "data in a tuberculosis clinical trial. Clin Transl Sci.",
    "2026;19:e70694. doi:10.1111/cts.70694.",
    "Final estimates from Supplementary Table S2; structure from Methods",
    "Equation 2 (exponent signs corrected, see the vignette). Age centring",
    "from the per-protocol baseline table of the primary REMoxTB report:",
    "Gillespie SH et al. Four-month moxifloxacin-based regimens for",
    "drug-sensitive tuberculosis. N Engl J Med. 2014;371(17):1577-1587.",
    "doi:10.1056/NEJMoa1407426.",
    sep = " "
  )
  vignette <- "You_2026_moxifloxacin_tb_ttp"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "h (ttp, MGIT time to positivity); no drug concentration"
  )

  # No d/dt() states: ttp is the algebraic prediction (canonical output in
  # inst/references/compartment-names.md).
  # One eta shared by both slopes (alpha_bi and beta_bi); there is no single
  # lslope_bi typical value for it to pair with.
  paper_specific_etas <- c("etalslope_bi")

  covariateData <- list(
    LUNG_CAVITY = list(
      description = "Cavitation on the baseline chest X-ray (1 = present, 0 = absent)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (cavitation present, the most common level)",
      notes = "Table S2 reports the effect of the ABSENCE of a cavity, so it enters as (1 + e_nocavity_abi * (1 - LUNG_CAVITY)). Source column CAVI (TB-PACTS RE domain REORRES, Y/N); missing values were replaced by the most frequent category (present).",
      source_name = "CAVI"
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the most common level)",
      notes = "Female is the non-reference level of the PsN stepwise covariate model.",
      source_name = "SEX"
    ),
    SMOKE = list(
      description = "Smoker (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-smoker)",
      notes = "Source column SMOKE = TB-PACTS SU domain SUOCCUR (Y = 1, N = 0). Read as current smoking: the effect is labelled 'Smoking on', i.e. smoking is the non-reference level, and only current smoking (31%) is less common than not smoking in the REMoxTB per-protocol population (Gillespie 2014 Table 1: never 46%, past 24%, current 31%).",
      source_name = "SMOKE"
    ),
    RACE_ASIAN = list(
      description = "Asian race (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other race categories)",
      notes = "Investigator-reported; includes South and East Asian (Gillespie 2014 Table 1).",
      source_name = "RACE ('ASIAN')"
    ),
    RACE_MISSING = list(
      description = "Race not reported (1 = not reported or 'other', 0 = race recorded)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (race recorded)",
      notes = "Table S2 'Unknown race'. The supplementary Python code pools the TB-PACTS RACE values 'NOT REPORTED' and 'OTHER' into one level, so both map to 1 here. With RACE_ASIAN = RACE_MISSING = 0 the subject is Black or mixed race.",
      source_name = "RACE ('NOT REPORTED', 'OTHER')"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 31.3 years, the grouped-data median of the REMoxTB per-protocol population (Gillespie 2014 Table 1: 30.2% < 25 y, 31.1% 25-35 y, 38.6% > 35 y; linear interpolation inside the 25-35 y band). You 2026 does not report its own analysis-dataset median, which PsN uses as the centre.",
      source_name = "AGE"
    ),
    TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL = list(
      description = "REMoxTB isoniazid group: moxifloxacin replaced ethambutol (isoniazid + rifampicin + pyrazinamide + moxifloxacin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control arm, isoniazid + rifampicin + pyrazinamide + ethambutol)",
      notes = "Table S2 pools this arm with the ethambutol group ('isoniazid or ethambutol'); the model uses the pooled moxifloxacin values when either indicator is 1.",
      source_name = "Arm"
    ),
    TRT_MOXIFLOXACIN_FOR_ISONIAZID = list(
      description = "REMoxTB ethambutol group: moxifloxacin replaced isoniazid (moxifloxacin + rifampicin + pyrazinamide + ethambutol)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control arm, isoniazid + rifampicin + pyrazinamide + ethambutol)",
      notes = "Table S2 pools this arm with the isoniazid group ('isoniazid or ethambutol').",
      source_name = "Arm"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1821L,
    n_studies = 1L,
    age_range = "Adults (>= 18 years); per-protocol population 30% < 25 y, 31% 25-35 y, 39% > 35 y (Gillespie 2014 Table 1)",
    weight_range = "Per-protocol population 10% < 40 kg, 16% 40-45 kg, 40% > 45-55 kg, 32% > 55-75 kg, 2% > 75 kg (Gillespie 2014 Table 1)",
    sex_female_pct = 30,
    race_ethnicity = c(Black = 44, Asian = 31, `Mixed race` = 25, Other = 0.2),
    disease_state = "Newly diagnosed, previously untreated, smear-positive, rifampicin- and fluoroquinolone-susceptible pulmonary tuberculosis.",
    dose_range = "REMoxTB (NCT00864383): control 8 weeks isoniazid + rifampicin + pyrazinamide + ethambutol then 18 weeks isoniazid + rifampicin; isoniazid group moxifloxacin replacing ethambutol and ethambutol group moxifloxacin replacing isoniazid, each 17 weeks then 9 weeks placebo. Doses adjusted to body weight (Gillespie 2014 Methods); no exposure term in the model.",
    regions = "South Africa, India, Tanzania, Kenya, Thailand, Malaysia, Zambia, China and Mexico (randomisation countries, Gillespie 2014 Results)",
    notes = "You 2026 Results 'Clinical Data': 1821 patients in the analysis dataset (scenario A); TTP data from days 0-14 for this model. Demographics are those of the REMoxTB per-protocol population (N = 1548, Gillespie 2014 Table 1) because You 2026 prints no demographics table. Cavitation 71%, current smokers 31%, HIV-positive 7% in that population."
  )

  ini({
    # Supplementary Table S2. Time in days, TTP in hours. A_bi and B_bi
    # are on the log scale: each term of TTP is exp(A_bi) or exp(B_bi).
    a_bi <- 3.43; label("A_bi, log amplitude of the first TTP term (log h)") # Table S2 'A_bi 3.43 (RSE 2%)'
    lalpha_bi <- log(0.0419); label("alpha_bi, first slope, control arm (1/day)") # Table S2 'alpha_bi control 0.0419 (RSE 8.3%)'
    lalpha_bi_moxi <- log(0.0477); label("alpha_bi, first slope, moxifloxacin arms (1/day)") # Table S2 'alpha_bi isoniazid or ethambutol 0.0477 (RSE 6%)'
    b_bi <- 4.54; label("B_bi, log amplitude of the second TTP term, control arm (log h)") # Table S2 'B_bi control 4.54 (RSE 0.7%)'
    b_bi_moxi <- 4.56; label("B_bi, log amplitude of the second TTP term, moxifloxacin arms (log h)") # Table S2 'B_bi isoniazid or ethambutol 4.56 (RSE 0.6%)'
    lbeta_bi <- log(0.0645); label("beta_bi, second slope, control arm (1/day)") # Table S2 'beta_bi control 0.0645 (RSE 2.8%)'
    lbeta_bi_moxi <- log(0.0737); label("beta_bi, second slope, moxifloxacin arms (1/day)") # Table S2 'beta_bi isoniazid or ethambutol 0.0737 (RSE 2%)'

    # Covariate effects, PsN stepwise-covariate linear form:
    # P * (1 + theta) for a binary covariate, P * (1 + theta * (x - centre))
    # for a continuous one, applied to A_bi and B_bi on their log scale.
    e_nocavity_abi <- 0.124; label("Absence of lung cavitation on A_bi (fractional change)") # Table S2 'Absence of cavity on A_bi 0.124 (RSE 13.5%)'
    e_sexf_abi <- -0.00554; label("Female sex on A_bi (fractional change)") # Table S2 'Female sex on A_bi -0.00554 (RSE 30.5%)'
    e_smoke_abi <- 0.061; label("Smoking on A_bi (fractional change)") # Table S2 'Smoking on A_bi 0.061 (RSE 23.8%)'
    e_asian_bbi <- -0.0251; label("Asian race on B_bi (fractional change)") # Table S2 'Asian race on B_bi -0.0251 (RSE 17.6%)'
    e_racemissing_bbi <- -0.0144; label("Unreported race on B_bi (fractional change)") # Table S2 'Unknown race on B_bi -0.0144 (RSE 35.2%)'
    e_sexf_bbi <- 0.0223; label("Female sex on B_bi (fractional change)") # Table S2 'Female sex on B_bi 0.0223 (RSE 21.9%)'
    e_racemissing_betabi <- 0.963; label("Unreported race on beta_bi (fractional change)") # Table S2 'Unknown race on beta_bi 0.963 (RSE 37.1%)'
    e_age_betabi <- -0.00489; label("Age on beta_bi (fractional change per year from 31.3 years)") # Table S2 'Age on beta_bi -0.00489 (RSE 21.1%)'
    e_smoke_betabi <- -0.0867; label("Smoking on beta_bi (fractional change)") # Table S2 'Smoking on beta_bi v -0.0867 (RSE 28.8%)'

    # IIV, variances (Table S2 footnote: 'inter-individual variability on
    # variance scale'). One omega is reported for 'slope' with one shrinkage
    # value although the model has two slopes; it is encoded as one eta
    # shared by alpha_bi and beta_bi.
    etaa_bi ~ 0.0427 # Table S2 'omega2 baseline 0.0427 (RSE 10.2%, shrinkage 27%)'
    etalslope_bi ~ 0.785 # Table S2 'omega2 slope 0.785 (RSE 9.7%, shrinkage 47%)'

    propSd <- 0.3132; label("Proportional residual error (fraction)") # Table S2 'sigma2 0.0981 (RSE 2.8%, shrinkage 14%)', a variance: sqrt(0.0981)
  })

  model({
    # Pooled moxifloxacin indicator: Table S2 gives one set of values for
    # the isoniazid and ethambutol groups together.
    moxi <- TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL + TRT_MOXIFLOXACIN_FOR_ISONIAZID

    alpha_tv <- exp(lalpha_bi) * (1 - moxi) + exp(lalpha_bi_moxi) * moxi
    b_tv <- b_bi * (1 - moxi) + b_bi_moxi * moxi
    beta_tv <- exp(lbeta_bi) * (1 - moxi) + exp(lbeta_bi_moxi) * moxi

    abi_i <- a_bi *
      (1 + e_nocavity_abi * (1 - LUNG_CAVITY)) *
      (1 + e_sexf_abi * SEXF) *
      (1 + e_smoke_abi * SMOKE) *
      exp(etaa_bi)
    bbi_i <- b_tv *
      (1 + e_asian_bbi * RACE_ASIAN) *
      (1 + e_racemissing_bbi * RACE_MISSING) *
      (1 + e_sexf_bbi * SEXF)
    alpha_i <- alpha_tv * exp(etalslope_bi)
    beta_i <- beta_tv *
      (1 + e_racemissing_betabi * RACE_MISSING) *
      (1 + e_age_betabi * (AGE - 31.3)) *
      (1 + e_smoke_betabi * SMOKE) *
      exp(etalslope_bi)

    # Equation 2 with positive exponents: TTP rises from about 127 h at
    # day 0 (Figure 3a); time is days after the start of treatment.
    ttp <- exp(abi_i) * exp(alpha_i * time) + exp(bbi_i) * exp(beta_i * time)

    ttp ~ prop(propSd)
  })
}
