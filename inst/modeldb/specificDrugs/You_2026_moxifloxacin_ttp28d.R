You_2026_moxifloxacin_ttp28d <- function() {
  description <- "Exponential-approach time-to-positivity (TTP) early-bactericidal-activity model for the first 28 days of treatment in adults with newly diagnosed smear-positive pulmonary tuberculosis in the REMoxTB trial (control isoniazid + rifampicin + pyrazinamide + ethambutol, or moxifloxacin replacing ethambutol or isoniazid). MGIT liquid-culture TTP (h) rises from an intercept towards intercept + asymptote at a first-order rate; the rate and the asymptote are estimated separately for each of the three arms. Covariates: lung cavitation and smoking on the intercept; age, lung cavitation, Asian or Black race, sex and body weight on the asymptote; baseline TTP, Black race and sex on the rate. Exponential IIV on intercept, rate and asymptote; proportional residual error. A $PRED regression on time with no exposure term, no differential equations and no dosing events."
  reference <- paste(
    "You H, Simonsson USH. Applying a pharmacometrics-enabled machine",
    "learning analysis to predict 2-month culture conversion using phase 2a",
    "data in a tuberculosis clinical trial. Clin Transl Sci.",
    "2026;19:e70694. doi:10.1111/cts.70694.",
    "Final estimates from Supplementary Table S3; structure from Methods",
    "Equation 3 (exponent sign corrected, see the vignette). Age and weight",
    "centring from the per-protocol baseline table of the primary REMoxTB",
    "report: Gillespie SH et al. Four-month moxifloxacin-based regimens for",
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

  covariateData <- list(
    LUNG_CAVITY = list(
      description = "Cavitation on the baseline chest X-ray (1 = present, 0 = absent)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (cavitation present, the most common level)",
      notes = "Table S3 reports the effects of the ABSENCE of cavities, so they enter as (1 + e * (1 - LUNG_CAVITY)). Source column CAVI (TB-PACTS RE domain REORRES, Y/N); missing values were replaced by the most frequent category (present).",
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
    RACE_BLACK = list(
      description = "Black race (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other race categories)",
      notes = "Investigator-reported. With RACE_ASIAN = RACE_BLACK = 0 the subject is of mixed race or another or unreported race.",
      source_name = "RACE ('BLACK')"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 31.3 years, the grouped-data median of the REMoxTB per-protocol population (Gillespie 2014 Table 1: 30.2% < 25 y, 31.1% 25-35 y, 38.6% > 35 y; linear interpolation inside the 25-35 y band). You 2026 does not report its own analysis-dataset median, which PsN uses as the centre.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 51.0 kg, the grouped-data median of the REMoxTB per-protocol population (Gillespie 2014 Table 1: 9.8% < 40 kg, 16.3% 40-45 kg, 40.1% > 45-55 kg; linear interpolation inside the 45-55 kg band). You 2026 does not report its own analysis-dataset median.",
      source_name = "WT"
    ),
    TTP_MGIT_BASE = list(
      description = "Baseline MGIT time to positivity",
      units = "day",
      type = "continuous",
      reference_category = NULL,
      notes = "Converted to hours in model() (x 24) and centred at 126.8 h, the median predicted baseline TTP of You 2026 Results; the paper does not report the median observed baseline TTP. The source value is the baseline-visit TTP, with the screening value substituted when the baseline value was missing or culture-negative (supplementary TTE_prepar_all.Rmd).",
      source_name = "baseline DVTTP"
    ),
    TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL = list(
      description = "REMoxTB isoniazid group: moxifloxacin replaced ethambutol (isoniazid + rifampicin + pyrazinamide + moxifloxacin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control arm, isoniazid + rifampicin + pyrazinamide + ethambutol)",
      notes = "Table S3 rows labelled 'isoniazid'.",
      source_name = "Arm"
    ),
    TRT_MOXIFLOXACIN_FOR_ISONIAZID = list(
      description = "REMoxTB ethambutol group: moxifloxacin replaced isoniazid (moxifloxacin + rifampicin + pyrazinamide + ethambutol)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control arm, isoniazid + rifampicin + pyrazinamide + ethambutol)",
      notes = "Table S3 rows labelled 'ethambutol'.",
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
    notes = "You 2026 Results 'Clinical Data': 1821 patients in the analysis dataset (scenario A); TTP data from days 0-28 for this model. Demographics are those of the REMoxTB per-protocol population (N = 1548, Gillespie 2014 Table 1) because You 2026 prints no demographics table. Cavitation 71%, current smokers 31%, HIV-positive 7% in that population."
  )

  ini({
    # Supplementary Table S3. Time in days, TTP in hours.
    la_exp <- log(117); label("A_exp, TTP intercept at day 0 (h)") # Table S3 'A_exp 117 (RSE 1.1%)'
    lalpha_exp <- log(0.0472); label("alpha_exp, rate of TTP rise, control arm (1/day)") # Table S3 'alpha_exp control 0.0472 (RSE 13.3%)'
    lalpha_exp_moxi_e <- log(0.0441); label("alpha_exp, rate of TTP rise, isoniazid group (moxifloxacin for ethambutol) (1/day)") # Table S3 'alpha_exp isoniazid 0.0441 (RSE 13.2%)'
    lalpha_exp_moxi_h <- log(0.0599); label("alpha_exp, rate of TTP rise, ethambutol group (moxifloxacin for isoniazid) (1/day)") # Table S3 'alpha_exp ethambutol 0.0599 (RSE 11.3%)'
    lb_exp <- log(341); label("B_exp, asymptotic TTP rise above the intercept, control arm (h)") # Table S3 'B_exp control 341 (RSE 8.8%)'
    lb_exp_moxi_e <- log(425); label("B_exp, asymptotic TTP rise, isoniazid group (moxifloxacin for ethambutol) (h)") # Table S3 'B_exp isoniazid 425 (RSE 9.1%)'
    lb_exp_moxi_h <- log(341); label("B_exp, asymptotic TTP rise, ethambutol group (moxifloxacin for isoniazid) (h)") # Table S3 'B_exp ethambutol 341 (RSE 7.2%)'

    # Covariate effects, PsN stepwise-covariate linear form:
    # P * (1 + theta) for a binary covariate, P * (1 + theta * (x - centre))
    # for a continuous one.
    e_nocavity_aexp <- 0.154; label("Absence of lung cavitation on A_exp (fractional change)") # Table S3 'Absence of cavities on A_exp 0.154 (RSE 16.9%)'
    e_smoke_aexp <- 0.0559; label("Smoking on A_exp (fractional change)") # Table S3 'Smoking on A_exp 0.0559 (RSE 33.1%)'
    e_age_bexp <- -0.00662; label("Age on B_exp (fractional change per year from 31.3 years)") # Table S3 'Age on B_exp -0.00662 (RSE 17.7%)'
    e_nocavity_bexp <- 0.129; label("Absence of lung cavitation on B_exp (fractional change)") # Table S3 'Absence of cavity on B_exp 0.129 (RSE 39.2%)'
    e_asian_bexp <- -0.121; label("Asian race on B_exp (fractional change)") # Table S3 'Asian race on B_exp -0.121 (RSE 26.5%)'
    e_sexf_bexp <- 0.401; label("Female sex on B_exp (fractional change)") # Table S3 'Female sex on B_exp 0.401 (RSE 34.4%)'
    e_wt_bexp <- 0.00793; label("Body weight on B_exp (fractional change per kg from 51.0 kg)") # Table S3 'Weight on B_exp 0.00793 (RSE 19.3%)'
    e_black_bexp <- 0.1; label("Black race on B_exp (fractional change)") # Table S3 'Black race on B_exp 0.1 (RSE 20.6%)'
    e_ttpbase_alphaexp <- 0.00388; label("Baseline TTP on alpha_exp (fractional change per h from 126.8 h)") # Table S3 'Baseline TTP on alpha_exp 0.00388 (RSE 18%)'
    e_black_alphaexp <- -0.184; label("Black race on alpha_exp (fractional change)") # Table S3 'Black race on alpha_exp -0.184 (RSE 29.2%)'
    e_sexf_alphaexp <- -0.336; label("Female sex on alpha_exp (fractional change)") # Table S3 'Female sex on alpha_exp -0.336 (RSE 28.6%)'

    # IIV, variances (Table S2/S3 footnotes: omega2 on the variance scale).
    etala_exp ~ 0.0788 # Table S3 'omega2 intercept 0.0788 (RSE 7.2%, shrinkage 18%)'
    etalalpha_exp ~ 0.34 # Table S3 'omega2 slope 0.34 (RSE 14.8%, shrinkage 44%)'
    etalb_exp ~ 0.227 # Table S3 'omega2 asymptote 0.227 (RSE 6.4%, shrinkage 31%)'

    propSd <- 0.3225; label("Proportional residual error (fraction)") # Table S3 'sigma2 0.104 (RSE 2.5%, shrinkage 12%)', a variance: sqrt(0.104)
  })

  model({
    ctrl <- 1 - TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL - TRT_MOXIFLOXACIN_FOR_ISONIAZID

    alpha_tv <- exp(lalpha_exp) * ctrl +
      exp(lalpha_exp_moxi_e) * TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL +
      exp(lalpha_exp_moxi_h) * TRT_MOXIFLOXACIN_FOR_ISONIAZID
    b_tv <- exp(lb_exp) * ctrl +
      exp(lb_exp_moxi_e) * TRT_MOXIFLOXACIN_FOR_ETHAMBUTOL +
      exp(lb_exp_moxi_h) * TRT_MOXIFLOXACIN_FOR_ISONIAZID

    ttpbase_h <- TTP_MGIT_BASE * 24

    aexp_i <- exp(la_exp) *
      (1 + e_nocavity_aexp * (1 - LUNG_CAVITY)) *
      (1 + e_smoke_aexp * SMOKE) *
      exp(etala_exp)
    alpha_i <- alpha_tv *
      (1 + e_ttpbase_alphaexp * (ttpbase_h - 126.8)) *
      (1 + e_black_alphaexp * RACE_BLACK) *
      (1 + e_sexf_alphaexp * SEXF) *
      exp(etalalpha_exp)
    bexp_i <- b_tv *
      (1 + e_age_bexp * (AGE - 31.3)) *
      (1 + e_nocavity_bexp * (1 - LUNG_CAVITY)) *
      (1 + e_asian_bexp * RACE_ASIAN) *
      (1 + e_sexf_bexp * SEXF) *
      (1 + e_wt_bexp * (WT - 51.0)) *
      (1 + e_black_bexp * RACE_BLACK) *
      exp(etalb_exp)

    # Equation 3 with a negative exponent: TTP rises from A_exp towards
    # A_exp + B_exp (Figure 3b); time is days after the start of treatment.
    ttp <- aexp_i + bexp_i * (1 - exp(-alpha_i * time))

    ttp ~ prop(propSd)
  })
}
