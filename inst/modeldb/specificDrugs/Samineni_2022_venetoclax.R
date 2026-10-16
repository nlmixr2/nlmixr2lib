Samineni_2022_venetoclax <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "first-order elimination for oral venetoclax given with R-CHOP to adults",
    "with relapsed/refractory or previously untreated non-Hodgkin lymphoma,",
    "mostly diffuse large B-cell lymphoma (the phase 1b/2 CAVALLI study,",
    "GO27878). Every structural, covariate, random-effect and residual value",
    "is carried FIXED from the legacy venetoclax population PK model (505",
    "subjects with CLL/SLL, NHL and healthy volunteers); Samineni 2022",
    "estimated only the two multiplicative CAVALLI study effects, x1.63 on",
    "CL/F and x2.33 on V2/F. Relative bioavailability depends on the meal",
    "(fasting / low-fat reference / moderate-fat / high-fat / fed with the",
    "meal unspecified) and on a fixed power function of the dose (reference",
    "400 mg). CL/F is lowered by strong and moderate CYP3A inhibitors and by",
    "OATP1B3 inhibitors; V2/F is higher in patients than in healthy",
    "volunteers and lower in women. Log-normal IIV on CL/F, V2/F and F1;",
    "combined proportional and additive residual error. Eight companion",
    "exposure-response logistic models in the Samineni_2022_venetoclax_*",
    "family take the steady-state nominal-dose AUC derived from this model",
    "as an input column."
  )
  reference <- paste(
    "Samineni D, Huang W, Gibiansky L, Ding H, Zhang R, Li C, Sinha A,",
    "Rajwanshi R, Humphrey K, Bazeos A, Salem AH, Miles D.",
    "Population Pharmacokinetics and Exposure-Response Analyses for",
    "Venetoclax in Combination with R-CHOP in Relapsed/Refractory and",
    "Previously Untreated Patients with Diffuse Large B Cell Lymphoma.",
    "Adv Ther. 2022;39(1):598-618. doi:10.1007/s12325-021-01919-z.",
    "Parameter values from Table 2; the legacy model is Jones AK et al.,",
    "AAPS J 2016;18(5):1192-1202, doi:10.1208/s12248-016-9927-9.",
    sep = " "
  )
  vignette <- "Samineni_2022_venetoclax"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L (= ug/mL)"
  )

  covariateData <- list(
    DOSE = list(
      description = "Administered venetoclax dose (mg) of the current dose record.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE register entry: drives the fixed dose",
        "nonlinearity on relative bioavailability, F1 x (DOSE / 400)^-0.18",
        "(Table 2 'F1,DOSE nonlinearity (400 mg as reference)' = -0.18).",
        "Supply it on dose records and on observation records. CAVALLI doses",
        "were 200 (n = 7), 400 (3), 600 (8) and 800 mg (214) once daily",
        "(Table 1). Eq. 1 of the source computes the exposure metric with",
        "the protocol-assigned (nominal) dose in this term."
      ),
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed-vs-fasted indicator for the dose record: 1 = dosed with food, 0 = dosed fasting.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 with FED_LOWFAT = 1 (low-fat meal) is the bioavailability reference",
      notes = paste(
        "Source column FOOD (five levels: fasting, low fat, moderate fat, high",
        "fat, fed). FED = 0 selects the Table 2 fasting factor 0.34. With",
        "FED = 1, exactly one of FED_LOWFAT / FED_MEDFAT / FED_HIGHFAT selects",
        "that meal; FED = 1 with all three 0 selects the 'FOOD=fed' level",
        "(meal type not recorded, factor 1.23). CAVALLI patients took",
        "venetoclax with the protocol meal, which the source classifies as",
        "low fat (Discussion), so set FED = 1 and FED_LOWFAT = 1 for them."
      ),
      source_name = "FOOD"
    ),
    FED_LOWFAT = list(
      description = "Low-fat-meal indicator for the dose record: 1 = dosed after a low-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (low-fat meal; the bioavailability reference of the legacy model, F1 factor 1)",
      notes = paste(
        "Source FOOD = low fat. Table 2 prints no low-fat row because it is",
        "the reference level; the same lineage (Dave 2019, Table 2 'Low fat on",
        "F1' = 1) states the low-fat reference explicitly. The CAVALLI",
        "protocol meals (23% and 3.5% of calories from fat) were both",
        "classified as low fat."
      ),
      source_name = "FOOD"
    ),
    FED_MEDFAT = list(
      description = "Moderate-fat-meal indicator for the dose record: 1 = dosed after a moderate-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Source FOOD = moderate fat. F1 x 1.31 against the low-fat reference (Table 2 'F1,FOOD=moderate fat').",
      source_name = "FOOD"
    ),
    FED_HIGHFAT = list(
      description = "High-fat-meal indicator for the dose record: 1 = dosed after a high-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Source FOOD = high fat. F1 x 1.43 against the low-fat reference (Table 2 'F1,FOOD=high fat').",
      source_name = "FOOD"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Strong CYP3A inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no strong CYP3A inhibitor; weak inhibitors are pooled into the reference)",
      notes = paste(
        "Source column C3AHIB (Table 1 levels: not administered, weak but no",
        "moderate or strong, moderate but no strong, strong). Strong -> CL/F x",
        "0.18. Mutually exclusive with CONMED_CYP3A4_INH_MOD. Weak inhibitors",
        "carry no effect. 6 of 228 CAVALLI patients (2.6%) had a strong",
        "inhibitor."
      ),
      source_name = "C3AHIB"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Moderate (but no strong) CYP3A inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no moderate or strong CYP3A inhibitor)",
      notes = "Source C3AHIB level 'moderate but no strong'. CL/F x 0.84. 65 of 228 CAVALLI patients (28.5%).",
      source_name = "C3AHIB"
    ),
    CONMED_OATP1B_INH = list(
      description = "OATP1B3 inhibitor coadministered: 1 = yes.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no OATP1B3 inhibitor)",
      notes = paste(
        "Source column OATP3HIB ('OATP1B3 inhibitor flag'). CL/F x 0.85. The",
        "source names OATP1B3 only; the register entry covers OATP1B1 and",
        "OATP1B3 jointly because clinical inhibitors are not selective",
        "between them. 3 of 228 CAVALLI patients (1.3%)."
      ),
      source_name = "OATP3HIB"
    ),
    STUDY_CAVALLI = list(
      description = "CAVALLI (GO27878) study indicator: 1 = subject from CAVALLI, 0 = subject from the legacy-model studies.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the eight phase 1/2 studies of the legacy model)",
      notes = paste(
        "Table 2 'CL/F, GO27878' = 1.63 and 'V2/F, GO27878' = 2.33, the only",
        "parameters estimated in Samineni 2022. Set STUDY_CAVALLI = 1 to",
        "reproduce the CAVALLI exposures; STUDY_CAVALLI = 0 gives the legacy",
        "model. The source attributes the shift to the lower-fat protocol",
        "meals and to population differences, not to a specific covariate."
      ),
      source_name = "GO27878"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator: 1 = healthy volunteer, 0 = patient with CLL, SLL or NHL.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy volunteer is the V2/F reference of the legacy model)",
      notes = paste(
        "Table 2 'V2/F, CLL/SLL/NHL' = 1.71: V2/F is 1.71-fold higher in",
        "patients than in healthy volunteers, applied as 1.71^(1 - DIS_HEALTHY).",
        "Every CAVALLI patient has DIS_HEALTHY = 0."
      ),
      source_name = "CLL/SLL/NHL"
    ),
    SEXF = list(
      description = "Sex: 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 2 'V2/F,SEX=female' = 0.68. 103 of 232 CAVALLI patients (44.4%) were female.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 232L,
    n_studies = 1L,
    n_observations = "1100 quantifiable venetoclax concentrations used in the analysis",
    age_range = "18-85 years (median 64, mean 61.2)",
    weight_range = "46.0-182 kg (median 77.0, mean 77.8)",
    sex_female_pct = 44.4,
    race_ethnicity = c(
      White = 71.1,
      Asian = 2.2,
      Black = 1.7,
      `Native Hawaiian or Pacific Islander` = 1.3,
      Unknown = 23.7
    ),
    disease_state = paste(
      "Relapsed/refractory or previously untreated non-Hodgkin lymphoma:",
      "DLBCL 218 (94.0%), follicular lymphoma 10 (4.3%), transformed",
      "lymphoma 2 (0.9%), other 2 (0.9%). ECOG PS 0/1/2: 44.4/40.5/15.1%.",
      "IPI 0-1/2/3/4-5: 11.2/37.1/29.7/22.0%. Stage IV 65.9%."
    ),
    dose_range = paste(
      "Venetoclax 200-800 mg orally once daily (800 mg in 92.2%) on days 4-10",
      "of cycle 1 and days 1-10 of cycles 2-8 of 21-day cycles, with",
      "rituximab (8 cycles) and CHOP (6-8 cycles); 32 phase 1b patients",
      "received obinutuzumab-CHOP instead."
    ),
    regions = "North America 34.5%, Western Europe 43.1%, Eastern Europe 15.9%, other 6.5%",
    co_medication = "CYP3A inhibitors: none 2.2%, weak only 66.7%, moderate but no strong 28.5%, strong 2.6%; OATP1B3 inhibitor 1.3%",
    notes = paste(
      "Baseline demographics from Samineni 2022 Table 1. The PK analysis",
      "dataset also contained the phase 1b G-CHOP arm. The fixed legacy",
      "parameters come from 505 subjects with R/R CLL/SLL, R/R NHL and",
      "healthy volunteers in eight phase 1/2 studies."
    )
  )

  ini({
    # Every value below is the Samineni 2022 Table 2 'Legacy model estimate'
    # column. The legacy parameters were held fixed in Samineni 2022 (Results:
    # 'All model parameters were fixed, including the inter-individual and
    # residual variability'); only the two GO27878 study effects were
    # estimated. Time in day, dose in mg, concentration in mg/L.
    lcl <- fixed(log(447)); label("Apparent clearance CL/F (L/day)") # Table 2 'CL/F (L/day)' = 447
    lvc <- fixed(log(118)); label("Apparent central volume of distribution V2/F (L)") # Table 2 'V2/F (L)' = 118
    lq <- fixed(log(97.2)); label("Apparent intercompartmental clearance Q/F (L/day)") # Table 2 'Q/F (L/day)' = 97.2
    lvp <- fixed(log(119)); label("Apparent peripheral volume of distribution V3/F (L)") # Table 2 'V3/F (L)' = 119
    lka <- fixed(log(3.72)); label("First-order absorption rate constant ka (1/day)") # Table 2 'ka (1/day)' = 3.72
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 of the reference condition (unitless)") # low-fat meal, 400 mg: reference level, no Table 2 row

    # Covariate effects on CL/F (multiplicative factors)
    e_conmed_cyp3a4_inh_strong_cl <- fixed(0.18); label("Multiplicative factor on CL/F with a strong CYP3A inhibitor (unitless)") # Table 2 'CL/F, strong CYP3A inhibitor' = 0.18
    e_conmed_cyp3a4_inh_mod_cl <- fixed(0.84); label("Multiplicative factor on CL/F with a moderate CYP3A inhibitor (unitless)") # Table 2 'CL/F, moderate CYP3A inhibitor' = 0.84
    e_conmed_oatp1b_inh_cl <- fixed(0.85); label("Multiplicative factor on CL/F with an OATP1B3 inhibitor (unitless)") # Table 2 'CL/F,OATP1B3 inhibitor' = 0.85
    e_study_cavalli_cl <- 1.63; label("Multiplicative CAVALLI (GO27878) study effect on CL/F (unitless)") # Table 2 'CL/F, GO27878' = 1.63; Results 95% CI 1.49-1.76

    # Covariate effects on V2/F (multiplicative factors)
    e_dis_patient_vc <- fixed(1.71); label("Multiplicative factor on V2/F for CLL/SLL/NHL patients vs healthy volunteers (unitless)") # Table 2 'V2/F, CLL/SLL/NHL' = 1.71
    e_sexf_vc <- fixed(0.68); label("Multiplicative factor on V2/F for women (unitless)") # Table 2 'V2/F,SEX=female' = 0.68
    e_study_cavalli_vc <- 2.33; label("Multiplicative CAVALLI (GO27878) study effect on V2/F (unitless)") # Table 2 'V2/F,GO27878' = 2.33; Results 95% CI 2.17-2.54

    # Food and dose effects on F1 (reference: low-fat meal, 400 mg)
    e_fasted_fdepot <- fixed(0.34); label("Multiplicative factor on F1 when fasting vs a low-fat meal (unitless)") # Table 2 'F1,FOOD=fasting' = 0.34
    e_medfat_fdepot <- fixed(1.31); label("Multiplicative factor on F1 with a moderate-fat meal vs a low-fat meal (unitless)") # Table 2 'F1,FOOD=moderate fat' = 1.31
    e_highfat_fdepot <- fixed(1.43); label("Multiplicative factor on F1 with a high-fat meal vs a low-fat meal (unitless)") # Table 2 'F1,FOOD=high fat' = 1.43
    e_fed_fdepot <- fixed(1.23); label("Multiplicative factor on F1 when fed with the meal type unrecorded vs a low-fat meal (unitless)") # Table 2 'F1,FOOD=fed' = 1.23
    e_dose_fdepot <- fixed(-0.18); label("Power exponent of (DOSE / 400 mg) on F1 (unitless)") # Table 2 'F1,DOSE nonlinearity (400 mg as reference)' = -0.18

    # Inter-individual variances (log-normal); no covariances are reported
    etalcl ~ fixed(0.15) # Table 2 'omega2 CL/F' = 0.15 (CV = 39.1%, shrinkage 19.0%)
    etalvc ~ fixed(0.21) # Table 2 'omega2 V2/F' = 0.21 (CV = 45.3%, shrinkage 27.7%)
    etalfdepot ~ fixed(0.10) # Table 2 'omega2 F1' = 0.10 (CV = 31.2%, shrinkage 28.3%)

    # Residual error: Table 2 prints variances; SD = sqrt(variance)
    propSd <- fixed(sqrt(0.22)); label("Proportional residual SD (fraction)") # Table 2 'sigma2 prop' = 0.22
    addSd <- fixed(sqrt(3.07e-7)); label("Additive residual SD (mg/L)") # Table 2 'sigma2 add' = 3.07 x 10^-7
  })

  model({
    cl <- exp(lcl + etalcl) *
      e_conmed_cyp3a4_inh_strong_cl^CONMED_CYP3A4_INH_STRONG *
      e_conmed_cyp3a4_inh_mod_cl^CONMED_CYP3A4_INH_MOD *
      e_conmed_oatp1b_inh_cl^CONMED_OATP1B_INH *
      e_study_cavalli_cl^STUDY_CAVALLI
    vc <- exp(lvc + etalvc) *
      e_dis_patient_vc^(1 - DIS_HEALTHY) *
      e_sexf_vc^SEXF *
      e_study_cavalli_vc^STUDY_CAVALLI
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)

    # Food factor on F1: fasting 0.34; low fat 1 (reference); moderate fat
    # 1.31; high fat 1.43; fed with no meal subtype recorded 1.23.
    f_meal <- FED_LOWFAT + FED_MEDFAT * e_medfat_fdepot + FED_HIGHFAT * e_highfat_fdepot +
      (1 - FED_LOWFAT - FED_MEDFAT - FED_HIGHFAT) * e_fed_fdepot
    f_food <- (1 - FED) * e_fasted_fdepot + FED * f_meal
    fdepot <- exp(lfdepot + etalfdepot) * f_food * (DOSE / 400)^e_dose_fdepot

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
