Savic_2017_rifapentine <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral rifapentine and its active metabolite 25-O-desacetyl rifapentine in adults with smear-positive pulmonary tuberculosis (TB Trials Consortium Studies 29 and 29X; Savic 2017). All parent elimination forms the metabolite, which has one-compartment disposition. Relative bioavailability falls with dose (600 / 900 / 1200-1500 mg vs 450 mg), is lower fasting or with a low-fat meal than with a high-fat meal, higher in White and Asian than Black participants and lower with HIV infection; apparent clearance falls with age and is lower in women; the fraction metabolised is 1.34-fold higher above 600 mg. Correlated IIV on CL/F, metabolite CL/F and F; combined additive + proportional residual error for each analyte."
  reference <- paste(
    "Savic RM, Weiner M, MacKenzie WR, Engle M, Whitworth WC, Johnson JL,",
    "Nsubuga P, Nahid P, Nguyen NV, Peloquin CA, Dooley KE, Dorman SE;",
    "Tuberculosis Trials Consortium of the Centers for Disease Control and Prevention.",
    "Defining the optimal dose of rifapentine for pulmonary tuberculosis:",
    "Exposure-response relations from two phase II clinical trials.",
    "Clin Pharmacol Ther. 2017;102(2):321-331. doi:10.1002/cpt.634."
  )
  vignette <- "Savic_2017_rifapentine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(
      analyte = "rifapentine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "rifapentine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_desrpt = list(
      analyte = "25-O-desacetyl rifapentine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on rifapentine CL/F centred on the median age of 31 years: CL = CL/F * (1 + e_age_cl * (AGE - 31)) (Savic 2017 Table 2 footnote a).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Women have 0.81-fold the CL/F of men (Savic 2017 Table 2 'Sex effect on CL (fraction in female vs. male)').",
      source_name = "Sex"
    ),
    RACE_WHITE = list(
      description = "White race indicator (1 = White)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference race is Black. Participants of other or unreported race (4% in Table 1) are not assigned an effect by the paper; RACE_WHITE = RACE_ASIAN = 0 gives them the Black reference.",
      source_name = "Race"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference race is Black (Savic 2017 Table 2).",
      source_name = "Race"
    ),
    HIV_POS = list(
      description = "HIV infection indicator (1 = HIV-infected)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies relative bioavailability by 0.85 (Savic 2017 Table 2).",
      source_name = "HIV"
    ),
    FED_HIGHFAT = list(
      description = "1 = rifapentine taken with a high-fat meal (> 27 g fat), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Mutually exclusive with FED_LOWFAT. The paper's bioavailability reference is the high-fat meal; fasting is FED_HIGHFAT = FED_LOWFAT = 0 (Savic 2017 Supplementary Methods 'Pharmacokinetic sampling': high fat > 27 g, lower fat 1 to 27 g, fasting).",
      source_name = "Food"
    ),
    FED_LOWFAT = list(
      description = "1 = rifapentine taken with a lower-fat meal (1 to 27 g fat), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Mutually exclusive with FED_HIGHFAT. Multiplies bioavailability by 0.83 relative to a high-fat meal (Savic 2017 Table 2).",
      source_name = "Food"
    ),
    DOSE_RIFAPENTINE_MG = list(
      description = "Rifapentine dose of the administration, in mg",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Selects the dose-level relative bioavailability (<= 450 mg: 1; <= 600 mg: 0.91; <= 900 mg: 0.83; > 900 mg: 0.74) and the 1.34-fold higher fraction metabolised above 600 mg (Savic 2017 Table 2). The studies gave 450, 600, 900, 1200 and 1500 mg (Table 3); other doses fall into the next-higher dose band, which is an extrapolation.",
      source_name = "Dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 405L,
    n_studies = 2L,
    age_range = "adults; median 31 years (IQR reported as 29.0-33.0) in the rifapentine group",
    age_median = "31 years",
    weight_range = "median 55.0 kg (IQR reported as 53.7-56.7) in the rifapentine group",
    weight_median = "55.0 kg",
    sex_female_pct = 29.4,
    race_ethnicity = c(Black = 60, White = 19, Asian = 16, Other = 0.7, NotReported = 4),
    disease_state = "Smear-positive pulmonary tuberculosis, intensive-phase multidrug therapy (rifapentine with isoniazid, pyrazinamide and ethambutol); 9% HIV-infected",
    dose_range = "Rifapentine 10 mg/kg 5 days/week fasting (Study 29) or 10, 15 or 20 mg/kg 7 days/week with food (Study 29X), given as 150 mg tablets; absolute doses 450-1500 mg once daily for 8 weeks",
    regions = "Brazil, Hong Kong, Kenya, Peru, South Africa, Spain, Vietnam, Uganda, United States",
    notes = "405 rifapentine-treated participants of the 657 in the PK/PD analysis (Savic 2017 Table 1). 2203 rifapentine and desacetyl rifapentine concentrations; 383 participants had samples (sparse 1-3 samples or intensive 7 samples between weeks 2 and 8; Supplementary Methods and Results). The interquartile ranges in Table 1 look like confidence intervals of the median and are quoted as printed."
  )

  ini({
    # Rifapentine (parent); Savic 2017 Table 2. Values are apparent (per F)
    # and refer to a 450 mg dose taken with a high-fat meal.
    lcl <- log(1.86)
    label("Rifapentine apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' = 1.86 (RSE 5%)
    lvc <- log(12.77)
    label("Rifapentine apparent volume V/F (L)") # Table 2 'V/F (L)' = 12.77 (RSE 5%)
    lka <- log(0.07)
    label("Rifapentine absorption rate constant ka (1/h)") # Table 2 'ka (h-1)' = 0.07 (RSE 3%)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability of a 450 mg dose with a high-fat meal (reference)") # Table 2 'Bioavailability of 450-mg dose with high fat food (reference)' = 1

    # 25-O-desacetyl rifapentine (metabolite); Table 2. All parent elimination
    # forms metabolite (Supplementary Figure S1), so the metabolite parameters
    # are apparent values divided by the unidentifiable fraction metabolised.
    lcl_desrpt <- log(1.91)
    label("Desacetyl rifapentine apparent clearance CLm/Fm (L/h)") # Table 2 'CLm/Fm (L/h)' = 1.91 (RSE 6%)
    lvc_desrpt <- log(8.83)
    label("Desacetyl rifapentine apparent volume Vm/Fm (L)") # Table 2 'Vm/Fm (L)' = 8.83 (RSE 12%)

    # Covariate effects; Table 2
    e_dose_600mg_fdepot <- 0.91
    label("Relative bioavailability of a 600 mg dose vs 450 mg (fraction)") # Table 2 = 0.91 (RSE 6%)
    e_dose_900mg_fdepot <- 0.83
    label("Relative bioavailability of a 900 mg dose vs 450 mg (fraction)") # Table 2 = 0.83 (RSE 6%)
    e_dose_1200mg_fdepot <- 0.74
    label("Relative bioavailability of a 1200 or 1500 mg dose vs 450 mg (fraction)") # Table 2 '1200-mg dose' = 0.74 (RSE 8%); Results: 'combined 1,200 and 1,500 mg'
    e_fasted_fdepot <- 0.72
    label("Relative bioavailability fasting vs high-fat meal (fraction)") # Table 2 'Fasting effect (vs. high fat)' = 0.72 (RSE 7%)
    e_fed_lowfat_fdepot <- 0.83
    label("Relative bioavailability with a low-fat meal vs high-fat meal (fraction)") # Table 2 'Effect of low-fat food (vs. high-fat)' = 0.83 (RSE 20%)
    e_race_white_fdepot <- 1.19
    label("Relative bioavailability White vs Black race (fraction)") # Table 2 = 1.19 (RSE 36%)
    e_race_asian_fdepot <- 1.51
    label("Relative bioavailability Asian vs Black race (fraction)") # Table 2 = 1.51 (RSE 16%)
    e_hiv_pos_fdepot <- 0.85
    label("Relative bioavailability HIV-infected vs uninfected (fraction)") # Table 2 = 0.85 (RSE 44%)
    e_age_cl <- -0.00379
    label("Linear change in CL/F per year of age from 31 years (1/year)") # Table 2 'Age effect on CL' = -0.00379 (RSE 50%), footnote a
    e_sexf_cl <- 0.81
    label("CL/F in women relative to men (fraction)") # Table 2 'Sex effect on CL (fraction in female vs. male)' = 0.81 (RSE 36%)
    e_dose_gt600mg_fm <- 1.34
    label("Fraction metabolised above 600 mg relative to 600 mg or less (fraction)") # Table 2 'Dose effect (> 600 mg) on fraction metabolized, Fm' = 1.34 (RSE 4%)

    # IIV; Table 2 CV% converted with omega^2 = log(1 + CV^2):
    # CL/F 40% -> 0.148420, CLm/Fm 44% -> 0.176974, F 36% -> 0.121864.
    # Covariances = r * sd_i * sd_j with Table 2 correlations
    # CL-CLm 0.69, CL-F 0.55, CLm-F 0.45.
    etalcl + etalcl_desrpt + etalfdepot ~ c(
      0.148420,
      0.111828, 0.176974,
      0.0739683, 0.0660852, 0.121864
    )

    # Residual error; Table 2
    propSd <- 0.19
    label("Rifapentine proportional residual error (fraction)") # Table 2 'Proportional residual error, rifapentine (CV%)' = 19
    addSd <- 1.61
    label("Rifapentine additive residual error (ug/mL)") # Table 2 'Additive residual error, rifapentine' = 1.61 ug/mL
    propSd_desrpt <- 0.14
    label("Desacetyl rifapentine proportional residual error (fraction)") # Table 2 'Proportional residual error, metabolite (CV%)' = 14
    addSd_desrpt <- 1.33
    label("Desacetyl rifapentine additive residual error (ug/mL)") # Table 2 'Additive residual error, metabolite' = 1.33 ug/mL
  })
  model({
    # Dose-level relative bioavailability; the paper estimated one factor per
    # studied dose (450, 600, 900, and 1200 pooled with 1500 mg).
    f_dose <- 1
    if (DOSE_RIFAPENTINE_MG > 450) f_dose <- e_dose_600mg_fdepot
    if (DOSE_RIFAPENTINE_MG > 600) f_dose <- e_dose_900mg_fdepot
    if (DOSE_RIFAPENTINE_MG > 900) f_dose <- e_dose_1200mg_fdepot

    # Food: high-fat meal is the reference; fasting is neither meal indicator.
    f_food <- e_fed_lowfat_fdepot^FED_LOWFAT *
      e_fasted_fdepot^((1 - FED_HIGHFAT) * (1 - FED_LOWFAT))

    f_cov <- e_race_white_fdepot^RACE_WHITE *
      e_race_asian_fdepot^RACE_ASIAN *
      e_hiv_pos_fdepot^HIV_POS

    fm <- 1
    if (DOSE_RIFAPENTINE_MG > 600) fm <- e_dose_gt600mg_fm

    cl <- exp(lcl + etalcl) * (1 + e_age_cl * (AGE - 31)) * e_sexf_cl^SEXF
    vc <- exp(lvc)
    ka <- exp(lka)
    cl_desrpt <- exp(lcl_desrpt + etalcl_desrpt)
    vc_desrpt <- exp(lvc_desrpt)

    kel <- cl / vc
    kel_desrpt <- cl_desrpt / vc_desrpt

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_desrpt) <- fm * kel * central - kel_desrpt * central_desrpt

    f(depot) <- exp(lfdepot + etalfdepot) * f_dose * f_food * f_cov

    # Dose in mg and volumes in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc_desrpt <- central_desrpt / vc_desrpt

    Cc ~ add(addSd) + prop(propSd)
    Cc_desrpt ~ add(addSd_desrpt) + prop(propSd_desrpt)
  })
}
