Mehta_2020_umeclidinium <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled umeclidinium in adults with COPD receiving FF/UMEC/VI single-inhaler triple therapy, FF/VI + UMEC or UMEC/VI, with body weight, age and current smoking on apparent inhaled clearance and body weight on apparent central volume"
  reference <- "Mehta R, Farrell C, Hayes S, Birk R, Okour M, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium Bromide/Vilanterol in Patients with Chronic Obstructive Pulmonary Disease. Clin Pharmacokinet. 2020;59(1):67-79. doi:10.1007/s40262-019-00794-w"
  vignette <- "Mehta_2020_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses (umeclidinium, 62.5 ug nominal) are entered in ug and
  # volumes are in L, so `Cc` is in ug/L == ng/mL. Mehta 2020 reports
  # concentrations and exposures in pg/mL and pg*h/mL (assay LLOQ 10 pg/mL,
  # Sect. 2.2); multiply `Cc` by 1000 to compare against the published values.

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F and V2/F normalised to 70 kg (Sect. 3.2 'typical value of UMEC CL/F was 149 L/h for a non-smoking subject with COPD aged 60 years and weighing 70 kg'). UMEC dataset median 71.8 kg, range 35.4-154 kg (Table 2).",
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 60 years (Sect. 3.2). UMEC dataset median 66 years, range 41-88 years (Table 2).",
      source_name = "age"
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-smoker (former smoker; Table 7 labels the levels Former / Current)",
      notes = "1 = current smoker (40% of the UMEC dataset, Table 2). CL/F multiplier 1.28 (Table 6).",
      source_name = "smoking status"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "umeclidinium", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 622,
    n_studies = 3,
    n_observations = 2589,
    age_range = "41-88 years",
    age_median = "66 years",
    weight_range = "35.4-154 kg",
    weight_median = "71.8 kg",
    bmi_median = "25.3 kg/m2 (range 14.4-49.2)",
    sex_female_pct = 27,
    race_ethnicity = c(
      `White/Caucasian/European` = 70,
      `Asian - East Asian` = 14,
      `Asian - Japanese` = 14,
      `African American/African` = 2,
      `White - Arabic/North African` = 0.3
    ),
    smoker_pct = 40,
    fev1_pct_predicted_median = "39.7% (range 12.6-79.6)",
    disease_state = "chronic obstructive pulmonary disease (symptomatic, at risk of exacerbation)",
    dose_range = "umeclidinium 62.5 ug once daily by oral inhalation as FF/UMEC/VI 100/62.5/25 ug (single inhaler), FF/VI 100/25 ug + UMEC 62.5 ug, or UMEC/VI 62.5/25 ug",
    regions = "multinational; all East Asian and Japanese subjects were resident in China, Japan or Korea",
    notes = "Pooled Phase III studies FULFIL (CTT116853), IMPACT (CTT116855) and 200812 (Table 1). Demographics from Table 2 (UMEC dataset). 13% of UMEC observations were below the 10 pg/mL LLOQ and were handled with the NONMEM M3 method; 33% of trough samples were BQL."
  )

  ini({
    lcl <- log(149); label("Apparent inhaled clearance CL/F, non-smoker, 60 years, 70 kg (L/h)") # Table 6 'CL/F (L/h)' 149 [138, 160], RSE 3.62%; Sect. 3.2 'typical value of UMEC CL/F was 149 L/h for a nonsmoking subject with COPD aged 60 years and weighing 70 kg'
    lvc <- log(1100); label("Apparent central volume V2/F, 70 kg (L)") # Table 6 'V2/F (L)' 1100 [1030, 1170], RSE 3.07%
    lq <- fixed(log(854)); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 6 'Q/F (L/h)' 854 fixed; Sect. 3.2 'Both Q/F and V3/F were fixed to previously estimated values'
    lvp <- fixed(log(16200)); label("Apparent peripheral volume V3/F (L)") # Table 6 'V3/F (L)' 16,200 fixed
    lka <- log(18.6); label("Absorption rate constant KA (1/h)") # Table 6 'KA (h-1)' 18.6 [16.2, 21.0], RSE 6.67%

    # Continuous covariates: TVPK = theta_pop * (COV/REF)^theta (Eq. 1). The
    # references are 70 kg and 60 years (Sect. 3.2 typical-value statement),
    # confirmed by the paper's own examples: (80/60)^-0.648 = 0.830 ('17%
    # lower'), (100/70)^0.580 = 1.230 ('23% higher'), (100/70)^0.797 = 1.329
    # ('33% higher'), (40/70)^0.797 = 0.640 ('36% lower').
    e_wt_cl <- 0.580; label("Power exponent on (WT/70) for CL/F (unitless)") # Table 6 'Body weight on CL/F' 0.580 [0.409, 0.751], RSE 15.0%
    e_age_cl <- -0.648; label("Power exponent on (AGE/60) for CL/F (unitless)") # Table 6 'Age on CL/F' -0.648 [-0.979, -0.317], RSE 26.1%
    e_wt_vc <- 0.797; label("Power exponent on (WT/70) for V2/F (unitless)") # Table 6 'Body weight on V2/F' 0.797 [0.614, 0.980], RSE 25.5%
    # Categorical covariate: TVPK = theta_pop * theta^CAT (Eq. 2); Table 6
    # prints theta (the multiplier), stored here as log(theta).
    e_smoke_cl <- log(1.28); label("Log-scale effect of current smoking on CL/F (unitless)") # Table 6 'Smoking effect on CL/F' 1.28 [1.13, 1.45], RSE 11.7%; Sect. 3.2 'CL/F increases by 28% in a subject with COPD who smoked'

    # IIV: Table 6 'IIV, CV%' converted with omega^2 = log(1 + CV^2).
    # No variances or correlations are printed; etas are independent.
    etalcl ~ 0.1329 # Table 6 CL/F IIV 37.7 CV% -> log(1 + 0.377^2)
    etalvc ~ 0.2352 # Table 6 V2/F IIV 51.5 CV% -> log(1 + 0.515^2)
    etalq ~ 0.3699 # Table 6 Q/F IIV 66.9 CV% -> log(1 + 0.669^2)
    etalvp ~ 0.4986 # Table 6 V3/F IIV 80.4 CV% -> log(1 + 0.804^2)
    etalka ~ 0.3524 # Table 6 KA IIV 65.0 CV% -> log(1 + 0.650^2)

    # Residual error not reported in Mehta 2020; declared at fixed(0).
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    cl <- exp(lcl + e_wt_cl * log(WT / 70) + e_age_cl * log(AGE / 60) + e_smoke_cl * SMOKE + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / 70) + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
