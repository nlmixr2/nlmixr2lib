Mehta_2020_vilanterol <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled vilanterol in adults with COPD receiving FF/UMEC/VI single-inhaler triple therapy, FF/VI + UMEC, FF/VI or UMEC/VI, with body weight on apparent inhaled clearance and current smoking on apparent central volume"
  reference <- "Mehta R, Farrell C, Hayes S, Birk R, Okour M, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium Bromide/Vilanterol in Patients with Chronic Obstructive Pulmonary Disease. Clin Pharmacokinet. 2020;59(1):67-79. doi:10.1007/s40262-019-00794-w"
  vignette <- "Mehta_2020_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses (vilanterol, 25 ug nominal) are entered in ug and
  # volumes are in L, so `Cc` is in ug/L == ng/mL. Mehta 2020 reports
  # concentrations and exposures in pg/mL and pg*h/mL (assay LLOQ 10 pg/mL,
  # Sect. 2.2); multiply `Cc` by 1000 to compare against the published values.

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 70 kg (Sect. 3.3 'typical value of VI CL/F was 73.5 L/h for a subject with COPD weighting 70 kg'). VI dataset median 71.9 kg, range 35.4-154 kg (Table 2).",
      source_name = "weight"
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-smoker (former smoker; Table 10 labels the levels Former / Current)",
      notes = "1 = current smoker (39% of the VI dataset, Table 2). V2/F multiplier 1.46 (Table 9).",
      source_name = "smoking status"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "vilanterol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 817,
    n_studies = 3,
    n_observations = 3331,
    age_range = "41-88 years",
    age_median = "66 years",
    weight_range = "35.4-154 kg",
    weight_median = "71.9 kg",
    bmi_median = "25.3 kg/m2 (range 14.4-49.2)",
    sex_female_pct = 29,
    race_ethnicity = c(
      `White/Caucasian/European` = 65,
      `Asian - East Asian` = 18,
      `Asian - Japanese` = 13,
      `African American/African` = 3,
      `White - Arabic/North African` = 0.2
    ),
    smoker_pct = 39,
    fev1_pct_predicted_median = "39.7% (range 6.65-79.6)",
    disease_state = "chronic obstructive pulmonary disease (symptomatic, at risk of exacerbation)",
    dose_range = "vilanterol 25 ug once daily by oral inhalation as FF/UMEC/VI 100/62.5/25 ug (single inhaler), FF/VI 100/25 ug + UMEC 62.5 ug, FF/VI 100/25 ug or UMEC/VI 62.5/25 ug",
    regions = "multinational; all East Asian and Japanese subjects were resident in China, Japan or Korea",
    notes = "Pooled Phase III studies FULFIL (CTT116853), IMPACT (CTT116855) and 200812 (Table 1). Demographics from Table 2 (VI dataset). 21% of VI observations were below the 10 pg/mL LLOQ and were handled with the NONMEM M3 method; 62% of trough samples were BQL."
  )

  ini({
    lcl <- log(73.5); label("Apparent inhaled clearance CL/F, 70 kg (L/h)") # Table 9 'CL/F (L/h)' 73.5 [69.7, 77.3], RSE 3.86%; Sect. 3.3 'typical value of VI CL/F was 73.5 L/h for a subject with COPD weighting 70 kg'
    lvc <- log(352); label("Apparent central volume V2/F, non-smoker (L)") # Table 9 'V2/F (L)' 352 [333, 371], RSE 3.44%; Sect. 3.3 'typical value of VI V2/F was 352 L for a non-smoking subject'
    lq <- log(242); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 9 'Q/F (L/h)' 242 [230, 254], RSE 4.96%
    lvp <- log(2250); label("Apparent peripheral volume V3/F (L)") # Table 9 'V3/F (L)' 2250 [1670, 2830], RSE 22.8%
    lka <- fixed(log(19.6)); label("Absorption rate constant KA (1/h)") # Table 9 'KA (h-1)' 19.6 fixed; Sect. 3.3 'The absorption rate constant was fixed to a previously estimated value'

    # Continuous covariate: TVPK = theta_pop * (WT/70)^theta (Eq. 1). Checks:
    # (40/70)^0.444 = 0.780 ('22% lower'), (100/70)^0.444 = 1.172 ('17%
    # higher'), Sect. 3.3.
    e_wt_cl <- 0.444; label("Power exponent on (WT/70) for CL/F (unitless)") # Table 9 'Body weight on CL/F' 0.444 [0.281, 0.607], RSE 20.2%
    # Categorical covariate: TVPK = theta_pop * theta^CAT (Eq. 2).
    e_smoke_vc <- log(1.46); label("Log-scale effect of current smoking on V2/F (unitless)") # Table 9 'Smoking effect on V2/F' 1.46 [1.34, 1.59], RSE 11.9%; Sect. 3.3 'VI V2/F was 46% higher'

    # IIV: Table 9 'IIV, CV%' converted with omega^2 = log(1 + CV^2).
    # KA carries IIV although its THETA is fixed (Table 9 lists 41.4 CV%).
    etalcl ~ 0.0797 # Table 9 CL/F IIV 28.8 CV% -> log(1 + 0.288^2)
    etalvc ~ 0.1807 # Table 9 V2/F IIV 44.5 CV% -> log(1 + 0.445^2)
    etalq ~ 0.0292 # Table 9 Q/F IIV 17.2 CV% -> log(1 + 0.172^2)
    etalvp ~ 0.6801 # Table 9 V3/F IIV 98.7 CV% -> log(1 + 0.987^2)
    etalka ~ 0.1582 # Table 9 KA IIV 41.4 CV% -> log(1 + 0.414^2)

    # Residual error not reported in Mehta 2020; declared at fixed(0).
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    cl <- exp(lcl + e_wt_cl * log(WT / 70) + etalcl)
    vc <- exp(lvc + e_smoke_vc * SMOKE + etalvc)
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
